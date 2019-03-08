//below are the final quality cuts I (Andrew Ludwig) used for the a4 analysis

#include "FFTtools.h"
#define CUTS_LOADED

TCut isRealImpulsive = "mostImpulsivePeak(2).value > 0 && flags.isRF == 1 && mostImpulsivePeak(2).theta < 60 && mostImpulsivePeak(2).theta > -50";
TCut isMCImpulsive = "abs(FFTtools::wrap(mostImpulsivePeak(2).phi - mc.phi,360,0)) < 5 && abs(FFTtools::wrap(mostImpulsivePeak(2).theta - mc.theta,360,0)) < 2.5";

TCut blastCutH = "!flags.isPayloadBlast && flags.maxBottomToTopRatio[0] < 2.8 && flags.maxBottomToTopRatio[0] > 0.98";
TCut blastCutV = "flags.maxBottomToTopRatio[1] < 2.8 && flags.maxBottomToTopRatio[1] > 0.98";
TCut nSectorsCut = "flags.nSectorsWhereBottomExceedsTop < 28";
TCut blastSum = "flags.maxBottomToTopRatio[0] + flags.maxBottomToTopRatio[1] < 4.3";
TCut totalBlastCut = blastSum && blastCutH && blastCutV && nSectorsCut;

TCut coherenceCut = "mostImpulsivePeak(2).antennaPeakAverage - mostImpulsiveCoherentFiltered(2).peakHilbert < 60"; 
TCut coherenceCut2 = "(mostImpulsivePeak(2).antennaPeakAverage/mostImpulsiveCoherentFiltered(2).peakHilbert*(mostImpulsivePeak(2).antennaPeakAverage - mostImpulsiveCoherentFiltered(2).peakHilbert)) < 250"; //changed from 200 for wais
TCut coherenceCut3 = "(mostImpulsivePeak(2).antennaPeakAverage/mostImpulsiveCoherentFiltered(2).peakHilbert) < 4"; //changed from 5
TCut totalCoherenceCut = coherenceCut && coherenceCut2 && coherenceCut3;

TCut stepFn = "flags.isStepFunction == 0 || flags.isStepFunction == 16";
TCut glitch = "flags.hasGlitch";
TCut fullGlitch = stepFn && !glitch;

TCut diCut = "mostImpulsiveDeconvolvedFiltered(2).fracPowerWindowGradient() <=20 || (flags.topPower[0]/flags.topPower[2] + flags.middleOrBottomPower[0]/flags.middleOrBottomPower[2]) < 1 / (mostImpulsiveDeconvolvedFiltered(2).fracPowerWindowGradient() - 20) + .575";
TCut usableCiCut = "(flags.topPower[0]/flags.topPower[2] + flags.middleOrBottomPower[0]/flags.middleOrBottomPower[2]) > .33 && (mostImpulsiveCoherentFiltered(2).fracPowerWindowGradient()) > 219.";
TCut squareCuts = !usableCiCut && diCut;

TCut lowPowerRatio = "flags.middleOrBottomPower[0] / flags.topPower[0] < 7.";
TCut allPowerRatio = "flags.middleOrBottomPower[2] / flags.topPower[2] < 3.5";
TCut allPowerRatio2 = "flags.topPower[2] / flags.middleOrBottomPower[2] < 6.";
TCut highPowerCut = "flags.middleOrBottomPower[1]/flags.middleOrBottomPower[2] < .6 && flags.topPower[1]/flags.topPower[2] < .7";
TCut powerRatioCuts = lowPowerRatio && allPowerRatio && allPowerRatio2 && highPowerCut; 

TCut qCuts = totalCoherenceCut  && isRealImpulsive && fullGlitch && totalBlastCut && powerRatioCuts && squareCuts; //these were my final sample qcuts
TCut signal_cut= isMCImpulsive && isRealImpulsive && fullGlitch && totalBlastCut && totalCoherenceCut && powerRatioCuts && squareCuts;

