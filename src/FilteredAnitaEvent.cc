#include "FilteredAnitaEvent.h" 
#include "UsefulAnitaEvent.h"
#include "FilterStrategy.h"
#include "RawAnitaHeader.h"
#include "AnalysisWaveform.h"
#include "AnitaVersion.h" 
#include "AnitaGeomTool.h"
#include <algorithm>
#include "TCanvas.h" 
#include "TStyle.h" 


//ClassImp(FilteredAnitaEvent); 

FilteredAnitaEvent::FilteredAnitaEvent() 
{
  for (unsigned i = 0; i < 2 * NUM_SEAVEYS; i++ )
    {
    filteredGraphs[i] = NULL; 
    rawGraphs[i] = NULL; 
  }
}

FilteredAnitaEvent:: FilteredAnitaEvent(const UsefulAnitaEvent * event, FilterStrategy * strategy, const Adu5Pat * rawpat, const RawAnitaHeader * header, bool save_stages ) 
  : useful(event), 
    strategy(strategy), 
    pat(rawpat), 
    header(header),
    keep_all_stages(save_stages)

{

#ifdef MULTIVERSION_ANITA_ENABLED  
  anitaVersion = AnitaVersion::getVersionFromUnixTime(header->realTime); 
#else 
  anitaVersion = 0; 
#endif 

  if (anitaVersion <= 0) anitaVersion = AnitaVersion::get(); 

  AnitaGeomTool * geom = AnitaGeomTool::Instance(); 
  // Initialize the filtered graphs with the raw graphs from Raw Anita Event 
  
  for (int pol = AnitaPol::kHorizontal; pol <= AnitaPol::kVertical; pol++)
  {
    for (unsigned ant = 0; ant < NUM_SEAVEYS; ant++) 
    {
      int k = pol * NUM_SEAVEYS  + ant; 
      int i = geom->getChanIndexFromAntPol(ant, (AnitaPol::AnitaPol_t) pol); 
      filteredGraphs[k] = new AnalysisWaveform (useful->fNumPoints[i], useful->fTimes[i], useful->fVolts[i]); 
      rawGraphs[k] = new AnalysisWaveform (useful->fNumPoints[i], useful->fTimes[i], useful->fVolts[i]); 
      filteredGraphs[k]->forceEvenSize(260); // do this for correlations 
//      rawGraphs[k]->forceEvenSize(260); // do this for correlations 
      filteredGraphsByAntPol[pol][ant] = filteredGraphs[k];  
      rawGraphsByAntPol[pol][ant] = rawGraphs[k];  
    }
  }


  //tell the strategy to process this
  strategy->process(this); 
}



const AnalysisWaveform * FilteredAnitaEvent::getFilteredGraphAtStage(UInt_t ant, AnitaPol::AnitaPol_t pol, UInt_t stage) const
{
  if (!keep_all_stages)
  {
    fprintf(stderr,"You didn't ask to save all the stages!\n"); 
    return 0; 
  }

  if (stage >= all_stages[pol][ant].size()) return filteredGraphsByAntPol[pol][ant]; 
  else return all_stages[pol][ant][stage]; 
}


void FilteredAnitaEvent::saveStage(int nreserve)
{

  for (int pol = 0; pol < 2; pol++)
  {
    for (int ant = 0; ant < NUM_SEAVEYS; ant++)
    {
      all_stages[pol][ant].reserve(nreserve); 
      all_stages[pol][ant].push_back(new AnalysisWaveform(*filteredGraphsByAntPol[pol][ant])); 
    }
  }
}

FilteredAnitaEvent::~FilteredAnitaEvent() 
{
  for (unsigned pol = 0; pol < 2; pol++)
  {
    for (unsigned ant = 0; ant <  NUM_SEAVEYS; ant++ )
    {
      delete filteredGraphsByAntPol[pol][ant]; 
      delete rawGraphsByAntPol[pol][ant];  

      for (unsigned j = 0; j < all_stages[pol][ant].size() ; j++)
      {
        delete all_stages[pol][ant][j]; 
      }
    }
  }
}


double FilteredAnitaEvent::getAveragePower(AnitaPol::AnitaPol_t pol, AnitaRing::AnitaRing_t ring, bool filtered) const
{

  AnitaPol::AnitaPol_t pol_start =pol == AnitaPol::kVertical ? AnitaPol::kVertical : AnitaPol::kHorizontal;  
  AnitaPol::AnitaPol_t pol_end = pol == AnitaPol::kHorizontal ? AnitaPol::kHorizontal : AnitaPol::kVertical; 

  double sum = 0; 

  int n = NUM_SEAVEYS * (int(pol_end) - int(pol_start) +1); 
  if (ring != AnitaRing::kNotARing) n/=AnitaRing::kNotARing; 

  for (int pol = (int) pol_start; pol <=(int) pol_end; pol++ ) 
  {
    for (int ant = 0; ant < NUM_SEAVEYS; ant++) 
    {
      if (ring != AnitaRing::kNotARing && AnitaGeomTool::getRingFromAnt(ant) != ring) continue; 
      sum += ( filtered ?  filteredGraphsByAntPol[pol][ant]->even() : rawGraphsByAntPol[pol][ant]->uneven())->getSumV2(); 
    }
  }

  return sum/n; 

}


double FilteredAnitaEvent::getMedianPower(AnitaPol::AnitaPol_t pol, AnitaRing::AnitaRing_t ring, bool filtered) const
{

  AnitaPol::AnitaPol_t pol_start =pol == AnitaPol::kVertical ? AnitaPol::kVertical : AnitaPol::kHorizontal;  
  AnitaPol::AnitaPol_t pol_end = pol == AnitaPol::kHorizontal ? AnitaPol::kHorizontal : AnitaPol::kVertical; 

  int n = NUM_SEAVEYS * (int(pol_end) - int(pol_start) +1); 

  if (ring != AnitaRing::kNotARing) n/=AnitaRing::kNotARing; 

  double vals[n]; 


  int i = 0; 
  for (int pol = (int) pol_start; pol <=(int) pol_end; pol++ ) 
  {
    for (int ant = 0; ant < NUM_SEAVEYS; ant++) 
    {
      if (ring != AnitaRing::kNotARing && AnitaGeomTool::getRingFromAnt(ant) != ring) continue; 
      vals[i++] = ( filtered ?  filteredGraphsByAntPol[pol][ant]->even() : rawGraphsByAntPol[pol][ant]->uneven())->getSumV2(); 
    }
  }

  std::nth_element(vals, vals + n/2, vals +n); 
  return vals[n/2];
}


void FilteredAnitaEvent::getMedianSpectrum(TGraph * target, AnitaPol::AnitaPol_t pol, double frac) const
{
  target->Set(131); 

  AnitaPol::AnitaPol_t pol_start =pol == AnitaPol::kVertical ? AnitaPol::kVertical : AnitaPol::kHorizontal;  
  AnitaPol::AnitaPol_t pol_end = pol == AnitaPol::kHorizontal ? AnitaPol::kHorizontal : AnitaPol::kVertical; 

//  memset(target->GetY(),0, target->GetN() * sizeof(double)); 
  memcpy(target->GetX(),filteredGraphsByAntPol[pol_start][0]->power()->GetX(), target->GetN() * sizeof(double)); 

  int n = NUM_SEAVEYS * (int(pol_end) - int(pol_start) +1); 


  double vals[n]; 
  int i = 0; 

  //Can paralellize this loop if it's helpful 
  for (int j = 0; j < 131; j++) 
  {
    i = 0; 
    for (int pol = (int) pol_start; pol <=(int) pol_end; pol++ ) 
    {
      for (int ant = 0; ant < NUM_SEAVEYS; ant++) 
      {
        //TODO: Is it faster just to put things into a set? 
        vals[i++] = filteredGraphsByAntPol[pol][ant]->powerdB()->GetY()[j];  
      }
    }
    std::nth_element(vals, vals + int(n*frac), vals +n); 
    target->GetY()[j] = vals[int(n*frac)]; 
  }
}


void FilteredAnitaEvent::getAverageSpectrum(TGraph * target, AnitaPol::AnitaPol_t pol) const
{
  target->Set(131); 

  AnitaPol::AnitaPol_t pol_start =pol == AnitaPol::kVertical ? AnitaPol::kVertical : AnitaPol::kHorizontal;  
  AnitaPol::AnitaPol_t pol_end = pol == AnitaPol::kHorizontal ? AnitaPol::kHorizontal : AnitaPol::kVertical; 

  memset(target->GetY(),0, target->GetN() * sizeof(double)); 
  memcpy(target->GetX(),filteredGraphsByAntPol[pol_start][0]->power()->GetX(), target->GetN() * sizeof(double)); 

  int n = NUM_SEAVEYS * (int(pol_end) - int(pol_start) +1); 

  for (int pol = (int) pol_start; pol <=(int) pol_end; pol++ ) 
  {
    for (int ant = 0; ant < NUM_SEAVEYS; ant++) 
    {
      for (int j = 0; j < 131; j++) 
      {
        target->GetY()[j] += filteredGraphsByAntPol[pol][ant]->power()->GetY()[j]/n;  
      }
    }
  }

  for (int j = 0; j < 131; j++) 
  {
    target->GetY()[j] = 10 * TMath::Log10(target->GetY()[j]); 
  }


}

void FilteredAnitaEvent::getMinMaxRatio(AnitaPol::AnitaPol_t pol, double * max_ratio, double * min_ratio, int* max_sector, int* min_sector, AnitaRing::AnitaRing_t ring1 , AnitaRing::AnitaRing_t ring2, int nth )  const
{
  
  double max = 0; 
  double min = 0; 
  int imax = -1; 
  int imin = -1; 

  for (int i = 0; i < NUM_PHI; i++) 
  {
    int ant1 = AnitaGeomTool::getAntFromPhiRing(i, ring1); 
    int ant2 = AnitaGeomTool::getAntFromPhiRing(i, ring2); 

//    printf("%d %d %d\n", pol, ant1,ant2); 
//    printf("%p\n", rawGraphsByAntPol[pol][ant1]); 
//    printf("%p\n", rawGraphsByAntPol[pol][ant2]); 
    double peak1 = rawGraphsByAntPol[pol][ant1]->uneven()->pk2pk(nth,nth);
    double peak2 = rawGraphsByAntPol[pol][ant2]->uneven()->pk2pk(nth,nth); 

    double ratio = peak1/peak2; 

    if ( imax < 0 || ratio > max ) 
    {
      imax = i; 
      max = ratio;
    }
    if ( imax < 0 || ratio < min ) 
    {
      imin = i; 
      min = ratio;
    }
  }

  if (max_ratio) *max_ratio = max; 
  if (min_ratio) *min_ratio = min; 
  if (max_sector) *max_sector= imax; 
  if (min_sector) *min_sector= imin; 


}


void FilteredAnitaEvent::plotSummary(TCanvas * ch, TCanvas * cv) const
{

  if (!ch) ch = new TCanvas("FilteredAnitaEvent_hpol","HPol", 1000,1000); 
  if (!cv) cv = new TCanvas("FilteredAnitaEvent_vpol","VPol", 1000,1000); 

  ch->Clear(); 
  cv->Clear(); 

  ch->Divide(8,NUM_SEAVEYS/8); 
  cv->Divide(8,NUM_SEAVEYS/8); 


  for (int i =0; i < NUM_SEAVEYS; i++) 
  {
    ch->cd(i+1); 
    getRawGraph(i, AnitaPol::kHorizontal)->drawUneven("",1); 
    getFilteredGraph(i, AnitaPol::kHorizontal)->drawEven("lsame",2); 

    cv->cd(i+1); 
    getRawGraph(i, AnitaPol::kVertical)->drawUneven("",1); 
    getFilteredGraph(i, AnitaPol::kVertical)->drawEven("lsame",2); 
  }

}


