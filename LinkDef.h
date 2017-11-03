#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class AnitaEventReconstructor+;
#pragma link C++ class AnitaEventSummary+;
#pragma link C++ class AnitaEventSummary::PointingHypothesis+;
#pragma link C++ class AnitaEventSummary::WaveformInfo+;
#pragma link C++ class AnitaEventSummary::ChannelInfo+;
#pragma link C++ class AnitaEventSummary::EventFlags+;
#pragma link C++ class AnitaEventSummary::SourceHypothesis+;
#pragma link C++ class AnitaEventSummary::MCTruth+;
#pragma link C++ class AnitaEventSummary::PayloadLocation+;
#pragma link C++ class FilteredAnitaEvent+;
#pragma link C++ class FilterStrategy+;
#pragma link C++ class AnalysisWaveform+;
#pragma link C++ class AnalysisWaveform::PowerCalculationOptions;
#pragma link C++ class PrettyAnalysisWaveform+;
#pragma link C++ class TGraphAligned+;
#pragma link C++ class CorrelationSummaryAnita4+;

#pragma link C++ class AnitaTemplateSummary+;
#pragma link C++ class AnitaTemplateSummary::SingleTemplateResult+;
#pragma link C++ class AnitaTemplateMachine+;

#pragma link C++ class AnitaNoiseSummary+;
#pragma link C++ class NoiseMonitor+;

#pragma link C++ class FilterOperation+;
#pragma link C++ class UniformFilterOperation+;
#pragma link C++ class ConditionalFilterOperation+;
#pragma link C++ class SimplePassBandFilter+;
#pragma link C++ class SimpleNotchFilter+;
#pragma link C++ class ALFASincFilter+;
#pragma link C++ class ALFAButterworthFilter+;
#pragma link C++ class ALFALanczosFilter+;
#pragma link C++ class HybridFilter+;
#pragma link C++ class SumDifferenceFilter+;
#pragma link C++ class DigitalFilterOperation+;
#pragma link C++ class GeometricFilter+;
#pragma link C++ class GaussianTaper; 
#pragma link C++ class DeglitchFilter+; 

#pragma link C++ namespace AnitaResponse+;
#pragma link C++ namespace impulsivity+;

#pragma link C++ namespace polarimetry+;
#pragma link C++ class polarimetry::StokesAnalysis;

#pragma link C++ class AnitaResponse::DeconvolutionMethod+;
#pragma link C++ class AnitaResponse::NaiveDeconvolution+;
#pragma link C++ class AnitaResponse::BandLimitedDeconvolution+;
#pragma link C++ class AnitaResponse::AbstractResponse+;
#pragma link C++ class AnitaResponse::Response+;
#pragma link C++ class AnitaResponse::ResponseManager+;
#pragma link C++ class AnitaResponse::CompositeResponse+;

#pragma link C++ class AnitaResponse::DeconvolveFilter+;
#pragma link C++ class AnitaResponse::WienerDeconvolution+; 
#pragma link C++ class AnitaResponse::AllPassDeconvolution+; 
#pragma link C++ class AnitaResponse::ImpulseResponseXCorr+; 
#pragma link C++ class AnitaResponse::DeconvolutionMethod+;

#pragma link C++ class AnitaEventFaker+; 


#endif


