#include "BlindDataset.h"

#include <fstream>
#include "TFile.h"
#include "TTree.h"


const BlindDataset::strategy defaultStrat = BlindDataset::kInsertedVPolEvents | BlindDataset::kInsertedHPolEvents;


BlindDataset::BlindDataset(strategy theStrat, int run, bool decimated, WaveCalType::WaveCalType_t cal, int anita_version) : AnitaDataset(run, decimated, cal, anita_version){
  setStrategy(theStrat);
  loadedBlindTrees = false;
  zeroPointers();
  loadBlindTrees();
}


BlindDataset::BlindDataset(int run, bool decimated, WaveCalType::WaveCalType_t cal, int anita_version) : AnitaDataset(run, decimated, cal, anita_version){
  std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ": initialised BlindDataset without specifying blinding."
	    << " Implementing default blinding, " << getDescription(BlindDataset::kDefault) << std::endl;
  setStrategy(BlindDataset::kDefault);
  zeroPointers();
  loadBlindTrees();
}



void BlindDataset::zeroPointers(){
  loadedBlindTrees = false;
  for(int pol=0; pol < AnitaPol::kNotAPol; pol++){
    fBlindHeadTree[pol] = NULL;
    fBlindEventTree[pol] = NULL;
  }
  fBlindFile = NULL;
}



TString BlindDataset::getDescription(BlindDataset::strategy strat){

  TString description = "Current strategy: ";

  if(strat == kNoBlinding){
    description = "No blinding. ";
  }

  if(strat & kInsertedVPolEvents){
    description += "VPol events inserted. ";
  }

  if(strat & kInsertedHPolEvents){
    description += "HPol events inserted. ";
  }

  return description;
}



BlindDataset::strategy BlindDataset::setStrategy(BlindDataset::strategy newStrat){
  theStrat = newStrat;
  return theStrat;
}








BlindDataset::strategy BlindDataset::getStrategy(){
  return theStrat;
}










void BlindDataset::loadBlindTrees() {

  std::cout << __PRETTY_FUNCTION__ << std::endl;

  if(!loadedBlindTrees){

    // prepare to put ROOT's current directory pointer back to what it was before we opened the blinding file in read mode.
    TDirectory* originalDir = gDirectory;

    std::cout << __PRETTY_FUNCTION__ << ": here 1" << std::endl;

    // zero internal pointers so can check we find everything.
    fBlindFile = NULL;

    char calibDir[FILENAME_MAX];
    char fileName[FILENAME_MAX];
    char *calibEnv=getenv("ANITA_CALIB_DIR");
    if(!calibEnv) {
      char *utilEnv=getenv("ANITA_UTIL_INSTALL_DIR");
      if(!utilEnv){
	sprintf(calibDir,"calib");
      }
      else{
	sprintf(calibDir,"%s/share/anitaCalib",utilEnv);
      }
    }
    else {
      strncpy(calibDir,calibEnv,FILENAME_MAX);
    }

    std::cout << __PRETTY_FUNCTION__ << ": here 2" << std::endl;


    // these are the fake events, that will be inserted in place of some min bias events
    sprintf(fileName, "%s/insertedDataFile.root", calibDir);
    fBlindFile = TFile::Open(fileName);

    gDirectory->ls();

    if(fBlindFile){

      TString polPrefix[AnitaPol::kNotAPol];
      polPrefix[AnitaPol::kHorizontal] = "HPol";
      polPrefix[AnitaPol::kVertical] = "VPol";

      for(int pol=0; pol < AnitaPol::kNotAPol; pol++){

	TString headTreeName = polPrefix[pol] + "HeadTree";
	fBlindHeadTree[pol] = (TTree*) fBlindFile->Get(headTreeName);

	TString eventTreeName = polPrefix[pol] + "EventTree";
	fBlindEventTree[pol] = (TTree*) fBlindFile->Get(eventTreeName);

	// If you found the data then prepare for data reading
	if(fBlindHeadTree[pol] && fBlindEventTree[pol]){

	  // std::cout << "This is the header tree \t" << fBlindHeadTree[pol]->GetName() << std::endl;
	  // fBlindHeadTree[pol]->GetListOfBranches()->Print();

	  // std::cout << "\n\n\n\n\n" << std::endl;

	  // std::cout << "This is the event tree \t" << fBlindEventTree[pol]->GetName() << std::endl;
	  // fBlindEventTree[pol]->GetListOfBranches()->Print();

	  // std::cout << "\n\n\n\n\n" << std::endl;

	  fBlindHeadTree[pol]->SetBranchAddress("header", &fBlindHeader[pol]);
	  fBlindEventTree[pol]->SetBranchAddress("event", &fBlindEvent[pol]);
	}
	else{
	  // complain if you can't find the data
	  std::cerr << "Error in " << __PRETTY_FUNCTION__ << ": "
		    << "fBlindHeadTree[" << pol << "] = " << fBlindHeadTree[pol] << ", "
		    << "fBlindEventTree[" << pol << "] = " << fBlindEventTree[pol] << std::endl;
	}
      }
    }
    else{
      std::cerr << "Error in " << __PRETTY_FUNCTION__ << ": "
		<< "Unable to find " << fileName << " for inserted event blinding." << std::endl;
    }

    std::cout << __PRETTY_FUNCTION__ << ": here 3" << std::endl;

    // these are the min bias event numbers to be overwritten, with the entry in the fakeEventTree
    // that is used to overwrite the event
    sprintf(fileName,"%s/anita3OverwrittenEventInfo.txt",calibDir);
    std::ifstream overwrittenEventInfoFile(fileName);
    char firstLine[180];
    overwrittenEventInfoFile.getline(firstLine,179);
    UInt_t overwrittenEventNumber;
    Int_t fakeTreeEntry, pol;
    Int_t numEvents = 0;
    while(overwrittenEventInfoFile >> overwrittenEventNumber >> fakeTreeEntry >> pol){
      eventsToOverwrite.push_back(overwrittenEventNumber);
      fakeTreeEntries.push_back(fakeTreeEntry);
      AnitaPol::AnitaPol_t thePol = AnitaPol::AnitaPol_t(pol);
      polarityOfEventToInsert.push_back(thePol);
      std::cout << overwrittenEventNumber << "\t" << fakeTreeEntry << "\t" << thePol << std::endl;
      // std::cout << overwrittenEventInfo.at(overwrittenEventInfo.size()-1).first << "\t" << overwrittenEventInfo.at(overwrittenEventInfo.size()-1).second << std::endl;
      numEvents++;
    }
    if(numEvents==0){
      std::cerr << "Warning in " << __FILE__ << std::endl;
      std::cerr << "Unable to find overwrittenEventInfo" << std::endl;
    }

    // put ROOT's current directory pointer back to what it was before we opened the blinding file in read mode.
    gDirectory = originalDir;
  }
  loadedBlindTrees = true;
}






/**
 * Loop through list of events to overwrite for a given polarisation and return the fakeTreeEntry we need to overwrite
 *
 * @param pol is the polarity to consider blinding
 * @param eventNumber is the eventNumber, obviously
 *
 * @return -1 if we don't overwrite, the entry in the fakeTree otherwise
 */
Int_t BlindDataset::needToOverwriteEvent(AnitaPol::AnitaPol_t pol, UInt_t eventNumber){

  Int_t fakeTreeEntry = -1;
  for(UInt_t i=0; i < polarityOfEventToInsert.size(); i++){
    if(polarityOfEventToInsert.at(i)==pol && eventNumber == eventsToOverwrite.at(i)){
      fakeTreeEntry = fakeTreeEntries.at(i);
      break;
    }
  }
  return fakeTreeEntry;
}





void BlindDataset::overwriteHeader(RawAnitaHeader* header, AnitaPol::AnitaPol_t pol, Int_t fakeTreeEntry){

  fBlindHeadTree[pol]->GetEntry(fakeTreeEntry);

  // Retain some of the header data for camouflage
  UInt_t realTime = header->realTime;
  UInt_t triggerTimeNs = header->triggerTimeNs;
  UInt_t eventNumber = header->eventNumber;
  Int_t run = header->run;
  Int_t trigNum = header->trigNum;
  UInt_t turfId = header->turfEventId;

  (*header) = (*fBlindHeader[pol]);

  header->realTime = realTime;
  header->triggerTimeNs = triggerTimeNs;
  header->eventNumber = eventNumber;
  header->run = run;
  header->trigNum = trigNum;
  header->turfEventId = turfId;


}




void BlindDataset::overwriteEvent(UsefulAnitaEvent* useful, AnitaPol::AnitaPol_t pol, Int_t fakeTreeEntry){

  fBlindEventTree[pol]->GetEntry(fakeTreeEntry);

  UInt_t eventNumber = useful->eventNumber;
  // std::cout << std::endl << std::endl;
  // std::cout << eventNumber << "\t" << fBlindEvent[pol]->eventNumber << std::endl;
  // std::cout << useful << "\t" << pol << std::endl;
  // std::cout << (useful->chipIdFlag[0] & 0x3) << "\t" << (fBlindEvent[pol]->chipIdFlag[0]&0x3) << std::endl;

  UInt_t surfEventIds[NUM_SURF] = {0};
  UChar_t wrongLabs[NUM_SURF*NUM_CHAN] = {0};
  UChar_t rightLabs[NUM_SURF*NUM_CHAN] = {0};
  for(int surf=0; surf < NUM_SURF; surf++){
    surfEventIds[surf] = useful->surfEventId[surf];

    for(int chan=0; chan < NUM_CHAN; chan++){
      const int chanIndex = surf*NUM_CHAN + chan;
      wrongLabs[chanIndex] = UChar_t(fBlindEvent[pol]->getLabChip(chanIndex));
      rightLabs[chanIndex] = UChar_t(useful->getLabChip(chanIndex));
      // std::cout << chanIndex << "\t" << chipIdFlags[chanIndex] << "\t" << (chipIdFlags[chanIndex] & 0x3) << std::endl;
    }
  }

  (*useful) = (*fBlindEvent[pol]);



  useful->eventNumber = eventNumber;
  for(int surf=0; surf < NUM_SURF; surf++){
    useful->surfEventId[surf] = surfEventIds[surf];

    // here we manually set the bits in the chipId flag that correspond to the lab chip
    // this ensures that as you click through magic display the LABS will still go A->B->C->D->A...
    for(int chan=0; chan < NUM_CHAN; chan++){

      const int chanIndex = surf*NUM_CHAN + chan;
      // std::cerr << useful->chipIdFlag[chanIndex] << "\t";
      useful->chipIdFlag[chanIndex] -= wrongLabs[chanIndex];
      // std::cerr << useful->chipIdFlag[chanIndex] << "\t";
      useful->chipIdFlag[chanIndex] += rightLabs[chanIndex];
      // std::cerr << useful->chipIdFlag[chanIndex] << std::endl;


    }
  }

}



/**
 * Blinding implementation for the headers
 *
 * @param force_load forces the tree load (see AnitaDataset)
 *
 * @return a pointer to the RawAnitaHeader, that may or may not have been overwritten
 */
RawAnitaHeader* BlindDataset::header(bool force_load){

  // load unblinded header as you would with AnitaDataset
  RawAnitaHeader* header = AnitaDataset::header(force_load);


  // This is the blinding implementation for the header

  if(theStrat & kInsertedVPolEvents){
    Int_t fakeTreeEntry = needToOverwriteEvent(AnitaPol::kVertical, header->eventNumber);
    if(fakeTreeEntry > -1){
      overwriteHeader(header, AnitaPol::kVertical, fakeTreeEntry);
    }
  }


  if(theStrat & kInsertedHPolEvents){
    Int_t fakeTreeEntry = needToOverwriteEvent(AnitaPol::kHorizontal, header->eventNumber);
    if(fakeTreeEntry > -1){
      overwriteHeader(header, AnitaPol::kHorizontal, fakeTreeEntry);
    }
  }
  return header;
}




/**
 * Blinding implementation for the headers
 *
 * @param force_load forces the tree load (see AnitaDataset)
 *
 * @return a pointer to the RawAnitaHeader, that may or may not have been overwritten
 */
UsefulAnitaEvent* BlindDataset::useful(bool force_load){

  // load unblinded header as you would with AnitaDataset
  UsefulAnitaEvent* useful = AnitaDataset::useful(force_load);

  // This is the blinding implementation for the header

  if(theStrat & kInsertedVPolEvents){
    Int_t fakeTreeEntry = needToOverwriteEvent(AnitaPol::kVertical, useful->eventNumber);
    if(fakeTreeEntry > -1){
      overwriteEvent(useful, AnitaPol::kVertical, fakeTreeEntry);
    }
  }


  if(theStrat & kInsertedHPolEvents){
    Int_t fakeTreeEntry = needToOverwriteEvent(AnitaPol::kHorizontal, useful->eventNumber);
    if(fakeTreeEntry > -1){
      overwriteEvent(useful, AnitaPol::kHorizontal, fakeTreeEntry);
    }
  }
  return useful;
}
