#include "BandwidthMeasure.h" 
#include "AnalysisWaveform.h" 
#include "TGraphAligned.h" 
#include <algorithm> 


double bandwidth::bandwidthMeasure(const AnalysisWaveform * wf, TGraph* testGraph) 
{
	const TGraphAligned* powdb = wf->power();
	std::vector<double> powers;
	double normconst = 0;
	for(int i = 0; i < powdb->GetN(); i++)
	{
		if(powdb->GetX()[i] < .18) continue;
		powers.push_back(powdb->GetY()[i]);
		normconst += powdb->GetY()[i];
	}
	int N = powers.size();
	for(int i = 0; i < N; i++) powers[i] = powers[i]/normconst;

	std::sort(powers.begin(), powers.end(), [](double a, double b) { 
			return b < a;
			});

	TGraph* gTemp = new TGraph(N);
	gTemp->SetPoint(0, 0, powers[0]);
	for(int i = 1; i < N; i++)
	{
		gTemp->SetPoint(i, i, powers[i] + gTemp->GetY()[i-1]);
	}
	if(testGraph)
	{
		delete testGraph;
		testGraph = (TGraph*) gTemp->Clone();
	}
	double cdf = 0;
	for(int i = 0; i < N; i ++) cdf += gTemp->GetY()[i]/double(N);
	delete gTemp;
	powers.clear();
	cdf = fabs(cdf-1)*2;
	return cdf;
}



