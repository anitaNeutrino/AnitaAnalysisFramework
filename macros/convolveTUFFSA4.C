#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TGraph.h"
#include "AnalysisWaveform.h"
#include "FFTtools.h"
#include <stdio.h>
#include <sys/stat.h>

void convolveTUFF(int ant=1, char tmb='B', char config='B', const char* outdir="")
{
	mkdir(Form("data/responses/A4ImpulseTUFFs/%s", outdir), 0777);
	TString impnameH;
	TString impnameV;
	TString outNameH;
	TString outNameV;
	if(ant < 10)
	{
		impnameH = Form("%s/share/AnitaAnalysisFramework/responses/A4noNotches/0%d%cH.imp", getenv("ANITA_UTIL_INSTALL_DIR"), ant, tmb);
		impnameV = Form("%s/share/AnitaAnalysisFramework/responses/A4noNotches/0%d%cV.imp", getenv("ANITA_UTIL_INSTALL_DIR"), ant, tmb);
		outNameH = Form("data/responses/A4ImpulseTUFFs/%s/0%d%cH.imp", outdir, ant, tmb);
		outNameV = Form("data/responses/A4ImpulseTUFFs/%s/0%d%cV.imp", outdir, ant, tmb);
	} 
	else
	{
		impnameH = Form("%s/share/AnitaAnalysisFramework/responses/A4noNotches/%d%cH.imp", getenv("ANITA_UTIL_INSTALL_DIR"), ant, tmb);
		impnameV = Form("%s/share/AnitaAnalysisFramework/responses/A4noNotches/%d%cV.imp", getenv("ANITA_UTIL_INSTALL_DIR"), ant, tmb);
		outNameH = Form("data/responses/A4ImpulseTUFFs/%s/%d%cH.imp", outdir, ant, tmb);
		outNameV = Form("data/responses/A4ImpulseTUFFs/%s/%d%cV.imp", outdir, ant, tmb);
	}

	TGraph tdH(impnameH);
	TGraph tdV(impnameV);
		
	AnalysisWaveform awfH(tdH.GetN(), tdH.GetY(), tdH.GetX()[1] - tdH.GetX()[0], 0);
	(void) awfH.freq();
	double dfH = awfH.deltaF();
	int NfreqH = awfH.Nfreq();
	FFTWComplex* freqH = awfH.updateFreq();
	AnalysisWaveform awfV(tdV.GetN(), tdV.GetY(), tdV.GetX()[1] - tdV.GetX()[0], 0);
	(void) awfV.freq();
	FFTWComplex* freqV = awfV.updateFreq();

	TFile convfile(Form("%s/share/AnitaAnalysisFramework/tuffModels/config%c.root", getenv("ANITA_UTIL_INSTALL_DIR"), config));
	TGraph* gReal = (TGraph*) convfile.Get("gReal");
	TGraph* gImag = (TGraph*) convfile.Get("gImag");
	double phaseShift = TMath::ATan2(gImag->Eval(0), gReal->Eval(0));

	//convolve w tuff response
	for(int i = 0; i < NfreqH; i++)
	{
		//if(dfH*i > 1.5) continue;
		double reTemp = gReal->Eval(dfH*i);
		double imTemp = gImag->Eval(dfH*i);
		FFTWComplex temp(reTemp, imTemp);
		temp.setMagPhase(temp.getAbs(), temp.getPhase() - phaseShift);
		freqH[i] *= temp;
		freqV[i] *= temp;
	}
	
	TGraph* gOutH = new TGraph(awfH.even()->GetN(), awfH.even()->GetX(), awfH.even()->GetY());
	TGraph* gOutV = new TGraph(awfV.even()->GetN(), awfV.even()->GetX(), awfV.even()->GetY());
	
	ofstream outH(outNameH.Data());
	ofstream outV(outNameV.Data());

	for(int i = 0; i < gOutH->GetN(); i++)
	{
		if(.1*i < 10)
		{
			outH << to_string(double(.1*i)).substr(0,3) << " " << gOutH->GetY()[i] << "\n";
			outV << to_string(double(.1*i)).substr(0,3) << " " << gOutV->GetY()[i] << "\n";
			continue;
		}
		if(.1*i < 100)
		{
			outH << to_string(double(.1*i)).substr(0,4) << " " << gOutH->GetY()[i] << "\n";
			outV << to_string(double(.1*i)).substr(0,4) << " " << gOutV->GetY()[i] << "\n";
			continue;
		}
		outH << to_string(double(.1*i)).substr(0,5) << " " << gOutH->GetY()[i] << "\n";
		outV << to_string(double(.1*i)).substr(0,5) << " " << gOutV->GetY()[i] << "\n";
	}

	outH.close();
	outV.close();
}

void makeAverage(const char* indir = "", const char* outdir="")
{
	mkdir(Form("data/responses/A4ImpulseTUFFs/%s", outdir), 0777);
	TString inf = Form("data/responses/A4ImpulseTUFFs/%s/", indir);
  char* dir = gSystem->ExpandPathName(inf.Data());
  void* dirp = gSystem->OpenDirectory(dir);
  const char* entry;
  const char* fname[100];
  Int_t n = 0;
  TString str;
  while((entry = (char*) gSystem->GetDirEntry(dirp)))
  {
    str = entry;
    if(str.EndsWith(".imp")) fname[n++] = gSystem->ConcatFileName(dir, entry);
  }
  TGraph* gOut = new TGraph(fname[0]);
  for(int i = 1; i < n; i++)
  {
    TGraph gTemp(fname[i]);
    for(int j=0; j < gOut->GetN(); j++)
    {
      gOut->GetY()[j] += gTemp.GetY()[j];
    }
  }
  for(int j=0; j < gOut->GetN(); j++)
  {
    gOut->GetY()[j] /= n;
  }
	TString oName = Form("data/responses/A4ImpulseTUFFs/%s/%s.imp", outdir, indir);
	ofstream outf(oName.Data());
	for(int i = 0; i < gOut->GetN(); i++)
	{
		if(.1*i < 10)
		{
			outf << to_string(double(.1*i)).substr(0,3) << " " << gOut->GetY()[i] << "\n";
			continue;
		}
		if(.1*i < 100)
		{
			outf << to_string(double(.1*i)).substr(0,4) << " " << gOut->GetY()[i] << "\n";
			continue;
		}
		outf << to_string(double(.1*i)).substr(0,5) << " " << gOut->GetY()[i] << "\n";
	}
  outf.close();
}

void makePow(const char* indir = "", const char* outdir="")
{
	mkdir(Form("data/responses/TUFFs/%s", outdir), 0777);
	TGraph* inf = new TGraph(Form("data/responses/TUFFs/%s/%s.imp", outdir, indir));
  TGraph* outf = FFTtools::makePowerSpectrum(inf);
  outf->SetName("power");
  TFile* f = new TFile(Form("data/responses/TUFFs/%s/%s.root", outdir, indir), "RECREATE");
  f->cd();
  outf->Write();
  f->Close();
  delete inf;
  delete outf;
}

void convolveTUFFSA4()
{
	mkdir("data/responses/A4ImpulseTUFFs", 0777);
	for(int i = 1; i < 17; i++)
	{
		convolveTUFF(i, 'B', 'B', "notches_260_375_0");
		convolveTUFF(i, 'M', 'B', "notches_260_375_0");
		convolveTUFF(i, 'T', 'B', "notches_260_375_0");
	}
  makeAverage("notches_260_375_0", "averages");
  makePow("notches_260_375_0", "averages");
	for(int i = 1; i < 17; i++)
	{
		convolveTUFF(i, 'B', 'J', "notches_250_375_0");
		convolveTUFF(i, 'M', 'J', "notches_250_375_0");
		convolveTUFF(i, 'T', 'J', "notches_250_375_0");
	}
  makeAverage("notches_250_375_0", "averages");
  makePow("notches_250_375_0", "averages");
	for(int i = 1; i < 17; i++)
	{
		convolveTUFF(i, 'B', 'C', "notches_260_0_460");
		convolveTUFF(i, 'M', 'C', "notches_260_0_460");
		convolveTUFF(i, 'T', 'C', "notches_260_0_460");
	}
  makeAverage("notches_260_0_460", "averages");
  makePow("notches_260_0_460", "averages");
	for(int i = 1; i < 17; i++)
	{
		convolveTUFF(i, 'B', 'A', "notches_260_0_0");
		convolveTUFF(i, 'M', 'A', "notches_260_0_0");
		convolveTUFF(i, 'T', 'A', "notches_260_0_0");
	}
  makeAverage("notches_260_0_0", "averages");
  makePow("notches_260_0_0", "averages");
	for(int i = 1; i < 17; i++)
	{
		convolveTUFF(i, 'B', 'O', "notches_260_365_0");
		convolveTUFF(i, 'M', 'O', "notches_260_365_0");
		convolveTUFF(i, 'T', 'O', "notches_260_365_0");
	}
  makeAverage("notches_260_365_0", "averages");
  makePow("notches_260_365_0", "averages");
	for(int i = 1; i < 17; i++)
	{
		convolveTUFF(i, 'B', 'G', "notches_260_385_0");
		convolveTUFF(i, 'M', 'G', "notches_260_385_0");
		convolveTUFF(i, 'T', 'G', "notches_260_385_0");
	}
  makeAverage("notches_260_385_0", "averages");
  makePow("notches_260_385_0", "averages");
	for(int i = 1; i < 17; i++)
	{
		convolveTUFF(i, 'B', 'P', "notches_260_375_460");
		convolveTUFF(i, 'M', 'P', "notches_260_375_460");
		convolveTUFF(i, 'T', 'P', "notches_260_375_460");
	}
  makeAverage("notches_260_375_460", "averages");
  makePow("notches_260_375_460", "averages");
}
