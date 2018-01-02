#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TGraph.h"
//#include "AnalysisWaveform.h"
//#include "FFTtools.h"
#include <stdio.h>
#include <sys/stat.h>

void convolveTUFFS_trigger()
{
	TString impname;
	TString outName;
	TFile impfile(Form("/home/keith/tuff_trigger_responses/icemc/data/Anita3_ImpulseResponseTrigger.root"));
	TGraph* gTrig = (TGraph*) impfile.Get("gTrigPath");

	outName = Form("%s/share/AnitaAnalysisFramework/responses/TUFFs/trigconfigP.imp", getenv("ANITA_UTIL_INSTALL_DIR"));

// To -do on gTrig but dont know how to do it with a pointer

	AnalysisWaveform awf(gTrig->GetN(), gTrig->GetY(), gTrig->GetX()[1] - gTrig->GetX()[0], 0);
	(void) awf.freq();
	double df = awf.deltaF();
	int Nfreq = awf.Nfreq();
	FFTWComplex* freq = awf.updateFreq();
//        const char config = "A";
	TFile convfile(Form("%s/share/AnitaAnalysisFramework/tuffModels/configP.root", getenv("ANITA_UTIL_INSTALL_DIR")));
	TGraph* gReal = (TGraph*) convfile.Get("gReal");
	TGraph* gImag = (TGraph*) convfile.Get("gImag");
	double phaseShift = TMath::ATan2(gImag->Eval(0), gReal->Eval(0));
        cout << "testing" << endl;
	//convolve w tuff response
	for(int i = 0; i < Nfreq; i++)
	{
		if(df*i > 1.5) continue;
		double reTemp = gReal->Eval(df*i);
		double imTemp = gImag->Eval(df*i);
		FFTWComplex temp(reTemp,imTemp);
		temp.setMagPhase(temp.getAbs(), temp.getPhase() - phaseShift);
		freq[i] *= temp;
	}
	
	TGraph* gOut = new TGraph(awf.even()->GetN(), awf.even()->GetX(), awf.even()->GetY());
	
	ofstream out(outName.Data());
        cout << outName << endl;

	for(int i = 0; i < gOut->GetN(); i++)
	{
		if(.1*i < 10)
		{
			out << to_string(double(.1*i)).substr(0,3) << " " << gOut->GetY()[i] << "\n";
			continue;
		}
		if(.1*i < 100)
		{
			out << to_string(double(.1*i)).substr(0,4) << " " << gOut->GetY()[i] << "\n";
			continue;
		}
		out << to_string(double(.1*i)).substr(0,5) << " " << gOut->GetY()[i] << "\n";
	}

	out.close();
}


