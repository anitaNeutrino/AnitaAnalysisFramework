#ifndef ANALYSIS_WAVEFORM_HH
#define ANALYSIS_WAVEFORM_HH

class FFTWComplex; 
#include "TGraphAligned.h" 
#include "Math/Interpolator.h"
#include "FFTWindow.h" 

/** 
 * \brief
 * This class is intended to be the main storage vessel for ANITA waveforms. It is similar in principle to RFWaveform from FFTtools. 
 * \endbrief
 *
 * This class holds 3 coupled versions of the waveform: 
 *
 *   - uneven: An unevenly sampled time-domain waveform 
 *   - even: An evenly sampled time-domain waveform
 *   - freq An evenly sampled frequency-domain waveform 
 *
 *   Transforming between uneven and even is accomplished via interpolation (which is not strictly invertible) 
 *   Transforming between even and freq is accomplished via FT (which is more or less invertible, up to machine precision)
 *
 *   Because the transformation between even and uneven is not reversible, it usually only makes sense to access
 *   or modify uneven prior to modifying even or freq. If an access to uneven is attempted after modifying even or freq, 
 *
 *
 *   There are two types of accessors, const ones which retrieve const data: 
 *     - uneven() 
 *     - even() 
 *     - freq() 
 *
 *
 *  And ones which allow modification of the data and will then recompute the other two versions the next time they're requested: 
 *   - updateUneven(); 
 *   - updateEven(); 
 *   - updateFreq(); 
 *
 *  Those have two versions, one where you modify the return value and one where you replace the value. 
 *
 *  In order to make it seem magic, there is a lot of internal state that is complicated to reason about. Hopefully there are no bugs...
 *
 * 
 *   Cosmin Deaconu <cozzyd@kicp.uchicago.edu> 
 *
 */  

class AnalysisWaveform 
{


  public: 

    /** Enable (or disable) a bunch of debugging crap. Only possible with ANITA_ANALYSIS_DEBUG defined to avoid slow branches */
    static void enableDebug(bool enable);  

    /** Enable (or disable) conversion from even to uneven. The default is that it's not allowed, which means
     * that as soon as some transformation is done on the evenly-sampled waveform or in the frequency domain, the uneven waveform
     * is set equal to the even waveform. */ 
    static void allowEvenToUnevenConversion(bool allow); 

    /** The interpolation method used to transform in between evenly-sampled and unevenly-sampled waveforms */ 
    enum InterpolationType
    {
      AKIMA, 
      SPARSE_YEN,    
      REGULARIZED_SPARSE_YEN    
    }; 

    /** The interpolation method used in e.g. evalEven */ 
    enum EvenEvaluationType
    {
      EVAL_LINEAR, 
      EVAL_AKIMA
    }; 

    /** The default type used. May be changed  by user */ 
    static InterpolationType defaultInterpolationType; 

    /** Interpolation options, relevant only for certain interpolators */ 
    struct InterpolationOptions 
    {

      InterpolationOptions() 
        : max_distance(8), regularization_order(0), mu(1e-3) {; } 

      int max_distance;  
      int regularization_order; 
      double mu; 
    }; 


    /** default power calculation options, may be changed by user */ 
    class PowerCalculationOptions
    {
      public: 
      PowerCalculationOptions() 
        : method(FFT), window_size(64)  { type = &FFTtools::GAUSSIAN_WINDOW; } 

      enum method
      {
        FFT, 
        BARTLETT,
        WELCH
      } method;  


      int window_size; 
      const FFTtools::FFTWindowType *type; 
    }; 

    /** Default interpolation options, may be changed by user */ 
    static InterpolationOptions defaultInterpolationOptions; 

    /* Constructors */ 

    /** Constructor for unevenly sampled waveform 
     * @param Nt the number of samples
     * @param x the time values
     * @param y the voltage values
     * @param nominal_dt the mean sample rate, defaults to ANITA rate
     * @param type  the InterpolationType used, defaults to defaultInterpolationType (which is probalby what you want to change to influence e.g. FilteredAnitaEvent). 
     * @param opt  the InterpolationOptions used, defaults to defaultInterpolationOptions (which is probably what you want to change to influence e.g. FilteredAnitaEvent). 
     *
     * */ 
    AnalysisWaveform(int Nt, const double * x, const double * y, double nominal_dt = 1./2.6, InterpolationType type = defaultInterpolationType, InterpolationOptions * opt  = &defaultInterpolationOptions);

    /** Constructor for evenly sampled waveform. Uneven waveform is set to even. */ 
    AnalysisWaveform(int Nt, const double * y, double dt, double t0);  

    /** Constructor from frequency domain (and uneven is set to be the same) */ 
    AnalysisWaveform(int Nt, const FFTWComplex * f, double df, double t0);  

    /** empty, even constructor */
    AnalysisWaveform(int Nt = 260, double dt=1./2.6, double t0=0);

    /** Copy constructor. Will not blindly copy everything (like it won't bother copying anything that will have to be recalculated on its next access */ 
    AnalysisWaveform(const AnalysisWaveform & other);  


    AnalysisWaveform & operator=(const AnalysisWaveform & other); 


    /** Computes the (circular) correlation (in the frequency domain) of the two waveforms as a new waveform. Note that if you want to
     * correlate two traces, they should be padded first. This does not pad them for you, but will complain if they are not! To turn off the nagging, see below function. It is also assumed the two are of the same length.
     *
     * By default is no normalization done at all, the frequency values are simply multiplied appropriately and then everything is divided by 1/scale. 
     *
     * If window_normalize is non-zero, then the total power of the overlapping window of at least size window_normalize at each point is used to scale. 
     * Smaller window sizes will use the smallest window size power ot normalize. An overall factor of scale is still applied.
     * 
     *
     **/ 
    static AnalysisWaveform * correlation(const AnalysisWaveform * A, const AnalysisWaveform * B, int npadfreq = 0, double scale =1 , int window_normalize = 0); 

    static AnalysisWaveform * convolution(const AnalysisWaveform * A, const AnalysisWaveform * B, int npadfreq = 0, double scale =1 ); 

    /* Enables or disables (check) and nagging on attemptoin correlation
     * without zero padding. Useful if you actually want to do a circular
     * correlation or if the check is too expensive for your tastes.
     *
     * Note that this only affects the current thread. 
     **/
    static void setCorrelationNag(bool nag); 


    /** Compute the (circular) autocorrelation of this waveform.
     *
     * If you want non-circular correlation, you can pass 1 (or more) to npadtime. 
     * No scaling is applied, but you can set a scale by passing a scale parameter. 
     * User responsible for deleting returned value; 
     **/ 

    AnalysisWaveform * autoCorrelation(int npadtime =1, int npadfreq = 0, double scale = 1); 

    /** Checks if the even waveform is zeropadded (by comparing second half to zero) */ 
    bool checkIfPaddedInTime() const; 

    virtual ~AnalysisWaveform(); 

    /* Constant accessors. If you coerce the compiler into allowing modification, coupled waveforms won't be updated */ 

    /** Constant accessor for uneven waveform. Dont' coerce into non-const unless you know what you're doing */ 
    const TGraphAligned * uneven() const ;

    /** Constant accessor for even waveform. Dont' coerce into non-const unless you know what you're doing */ 
    const TGraphAligned * even()  const; 

    /** Constant accessor for frequency domain waveform. Dont' coerce into non-const unless you know what you're doing */ 
    const FFTWComplex * freq() const; 

    /** Get the power */ 
    const TGraphAligned * power() const; 

    /** Get the power in dB*/ 
    const TGraphAligned * powerdB() const; 

    /** Get the phase */ 
    const TGraphAligned * phase() const; 

    /** Get the group delay  */ 
    const TGraphAligned * groupDelay() const; 

    /** Get the hilbert envelope */ 
    const TGraphAligned * hilbertEnvelope() const; 

    /** Get the Hilbert Transform */ 
    const AnalysisWaveform * hilbertTransform() const; 

    /** Return the length of the frequency domain complex waveform */ 
    int Nfreq() const;

    /** Return the number of samples in the even waveform */ 
    int Neven() const { return even()->GetN(); } 

    /** Return the number of samples in the uneven waveform */ 
    int Nuneven() const { return uneven()->GetN(); } 

    /** Return the mean sampling period */ 
    double deltaT() const { return dt ; } 

    /** Return the spacing between frequencies in the frequency domain*/ 
    double deltaF() const { return df ; }

    /** Forces the even waveform to be a particular size. This is equivalent to 
     * updateEven()->Set() (with padding) with the exception that it may be deferred
     * if even waveform needs to be (re)calculated.*/
    void forceEvenSize(int size); 

    // drawers since drawing is non-const (and we don't care about silly things like axes for constness)
    
    /** Draw even */ 
    void drawEven(const char * opt = "", int color=-1) const; 

    /** Draw Hilbert envelope */ 
    void drawHilbertEnvelope(const char * opt = "", int color=-1) const; 

    /** Draw uneven */ 
    void drawUneven(const char * opt = "", int color=-1) const; 

    /** Draw power */
    void drawPower(const char * opt = "", int color=-1) const; 


    /** Draw power in dB*/
    void drawPowerdB(const char * opt = "", int color=-1) const; 

    /** Draw phase */
    void drawPhase(const char * opt = "", int color=-1) const; 

    /**Evaluate the even waveform at a point. */ 
    double evalEven(double t, EvenEvaluationType = EVAL_LINEAR) const; //TODO: add additional evaluation methods other than linear interpolation. indeed, would be best to eval multiple points at same time

    /** Evaluate the even waveform at the specified points. */ 
    void evalEven(int N, const double * __restrict t , double * __restrict v, EvenEvaluationType = EVAL_LINEAR) const; 


    /**  Update frequency graph by modifying return value */ 
    FFTWComplex * updateFreq(); 
    /**  Update frequency graph by replacing with argumetnts*/ 
    void updateFreq(int new_N, const FFTWComplex * new_freq, double new_df = 0) ;

    /** Update the even graph by modifying return value */ 
    TGraphAligned * updateEven(); 

    /** Update the even graph by replacing with argument*/
    void updateEven(const TGraph * replace_even); 

    /** Update the uneven graph by modifying return value */ 
    TGraphAligned * updateUneven(); 

    /** Update the uneven graph by replacing with argument*/
    void updateUneven(const TGraph * replace_uneven); 


    /** pad the even waveform (equivalent to upsampling the frequency). 
     *
     *
     * if where > 0, zeroes will be added at the end
     * if where < 0, zeroes will be added at the beginning
     * if where == 0, zeroes will be padded before and after
     *  
     * */
    void padEven(int factor, int where = 1); 

    /** pad the frequency (equivalent to upsampling the even values) */
    void padFreq(int factor); 


    /** pad the frequency with a specific number of points (equivalent to upsampling the even values) */
    void padFreqAdd(int npad);


    /** Change the options for calculating the power. This obviously will invalidate power() and powerdB() */
    void setPowerCalculationOptions(PowerCalculationOptions & opt); 


    /** change between vpol / hpol and lcp / rcp 
     *
     * */ 
    static void basisChange(AnalysisWaveform * __restrict hpol_or_lcp, AnalysisWaveform * __restrict vpol_or_rcp); 

    /** Change between a ,b -> (a + b)/2, (a - b)/2 
     * */ 
    static void sumDifference(AnalysisWaveform * __restrict a, AnalysisWaveform * __restrict b); 


    /** Zero mean even graph */ 
    void zeroMeanEven(); 

    void setTitle(const char * title); 
    /* Sets the line and marker color */ 
    void setColor(int c); 
    /* Sets the line width */ 
    void setWidth(int w);  

    /* Sets range user for x axis of frequency plots */ 
    void setFreqDisplayRange(double low, double high); 

    /* Sets range user for x axis of time plots */ 
    void setTimeDisplayRange(double low, double high); 

  private: 
    void calculateEvenFromUneven() const; 
    void calculateEvenFromFreq() const;  
    void calculateFreqFromEven() const;  
    void calculateUnevenFromEven() const;

    /* Storage vars */ 
    mutable TGraphAligned g_uneven; 
    mutable TGraphAligned g_hilbert_envelope; 
    mutable TGraphAligned g_even; 
    mutable TGraphAligned g_power; 
    mutable TGraphAligned g_power_db;
    mutable TGraphAligned g_phase;  
    mutable TGraphAligned g_group_delay; 


    mutable double dt;  
    mutable double df; 
    mutable int fft_len; 
    mutable FFTWComplex * fft; 

    mutable ROOT::Math::Interpolator theEvenAkimaInterpolator; 
    const ROOT::Math::Interpolator * evenAkimaInterpolator() const; 

    InterpolationType interpolation_type; 
    InterpolationOptions interpolation_options; 
    PowerCalculationOptions power_options; 

    /* State vars */ 
    mutable bool just_padded;  
    mutable bool must_update_uneven; 
    mutable bool must_update_freq; 
    mutable bool must_update_even; 
    mutable bool uneven_equals_even; 
    mutable bool power_dirty; 
    mutable bool power_db_dirty; 
    mutable bool phase_dirty; 
    mutable bool group_delay_dirty;
    mutable bool hilbert_envelope_dirty; 
    mutable bool hilbert_dirty; 
    mutable bool even_akima_interpolator_dirty; 
    mutable AnalysisWaveform * hilbert_transform;
    mutable int force_even_size; 
    unsigned uid; 
    void nameGraphs(); 
}; 



#endif 
