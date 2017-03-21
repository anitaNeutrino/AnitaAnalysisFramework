/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
 Separate the blinding implementation from AnitaEventCalibrator, in fact get it out of eventReaderRoot.
 I think is really is the natural place for it, since it is part of analysis, not a part of the data.
 To add a blinding strategy please update the strategy enum, and the description function in cc file.
 You will also need to add a function to manipulate the UsefulAnitaEvent you want to blind.
***********************************************************************************************************/

#ifndef BLIND_DATASET_H
#define BLIND_DATASET_H

#include <iostream>
#include "UsefulAnitaEvent.h"
#include "RawAnitaHeader.h"

#include "TString.h"
#include "AnitaDataset.h"

class BlindDataset : public AnitaDataset {

public:

  // All the blinding options.
  // The plan is to have each one represented by a bit so multiple strategies can be set at the same time.
  // Any output files should probably save the blinding strategy to the tree.
  enum strategy {
    kNoBlinding                          = 0x00,
    kInsertedVPolEvents                  = 0x01,
    kInsertedHPolEvents                  = 0x02,
    // kYetAnotherStrategy                  = 0x04,
    // kYetAnotherStill                     = 0x08


    kDefault= kInsertedVPolEvents|kInsertedHPolEvents
  };


  BlindDataset(strategy theStrat, int run, bool decimated, WaveCalType::WaveCalType_t cal = WaveCalType::kDefault, int anita_version = -1);
  BlindDataset(int run, bool decimated, WaveCalType::WaveCalType_t cal = WaveCalType::kDefault, int anita_version = -1);




  /**
   * Get a one line description of the blinding strategy
   * (please update this in BlindDataset.cc when you add a strategy.
   * @param strat is the strategy to describe.
   * @return the description
   */
  TString getDescription(strategy strat);


  /**
   * Set the current strategy (see BlindDataset.h) for strategy options
   * @param strat is the strategy to use. This can be a combination e.g. (kIntertedEvents | kAnotherStrategy)
   * @return the strategy that was set
   */
  strategy setStrategy(strategy strat);


  /**
   * Get the currently set strategy
   * @return the strategy that was set
   */
  strategy getStrategy();






  RawAnitaHeader* header(bool force_load=false); // add a few sneaky tricks to AnitaDataset::header(bool)
  UsefulAnitaEvent* useful(bool force_load=false); // add a few sneaky tricks to AnitaDataset::header(bool)



private:

  void zeroPointers();
  void loadBlindTrees();


  strategy theStrat; ///!< Currently selected blinding strategy
  bool loadedBlindTrees; ///!< Have we loaded the tree of events to insert?
  Int_t needToOverwriteEvent(AnitaPol::AnitaPol_t pol, UInt_t eventNumber);
  void overwriteHeader(RawAnitaHeader* header, AnitaPol::AnitaPol_t pol, Int_t fakeTreeEntry);
  void overwriteEvent(UsefulAnitaEvent* useful, AnitaPol::AnitaPol_t pol, Int_t fakeTreeEntry);

  // fake things
  TFile* fBlindFile; ///!< Pointer to file containing tree of UsefulAnitaEvents to insert
  TTree* fBlindEventTree[AnitaPol::kNotAPol]; ///!< Tree of UsefulAnitaEvents to insert
  TTree* fBlindHeadTree[AnitaPol::kNotAPol]; ///!< Tree of headers to insert

  UsefulAnitaEvent* fBlindEvent[AnitaPol::kNotAPol]; ///!< Pointer to fake UsefulAnitaEvent
  RawAnitaHeader* fBlindHeader[AnitaPol::kNotAPol]; ///!< Pointer to fake header

  // Here we have a set of three vectors that should all be the same length as the elements of each match up.
  // They are filled in loadBlindTrees
  std::vector<UInt_t> eventsToOverwrite;
  std::vector<AnitaPol::AnitaPol_t> polarityOfEventToInsert;
  std::vector<Int_t> fakeTreeEntries;

};





// define the bitwise or operator | to combine blinding strategies
inline BlindDataset::strategy operator|(BlindDataset::strategy strat1,  BlindDataset::strategy strat2){
  return static_cast<BlindDataset::strategy>(static_cast<UInt_t>(strat1) | static_cast<UInt_t>(strat2));
}
// define the bitwise and opterator & to decode the blinding strategies
inline BlindDataset::strategy operator&(BlindDataset::strategy strat1,  BlindDataset::strategy strat2){
  return static_cast<BlindDataset::strategy>(static_cast<UInt_t>(strat1) & static_cast<UInt_t>(strat2));
}

#endif //BLIND_DATASET_H
