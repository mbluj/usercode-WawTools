#ifndef WarsawAnalysis_HTTDataFormats_HTTEvent_h
#define WarsawAnalysis_HTTDataFormats_HTTEvent_h

#include "TLorentzVector.h"
#include "TVector3.h"

//
// class declaration
//

class HTTEvent {

 public:

  HTTEvent();

  ~HTTEvent();

  ///Reset data members to default values.
  void clear();

  
  float run_, lumi_, event_;
  int bosonId_, decModeMinus_, decModePlus_;

  int pfPVIndex_, pt2PVindex_;

  bool isRefit_;
  int nTracksInRefit_;
  
  TLorentzVector p4Sum_;
  TLorentzVector piMinus_, piPlus_;
  TLorentzVector tauMinus_, tauPlus_;
  TLorentzVector visTauMinus_, visTauPlus_;
  
  TVector3 genPV_, aodPV_, pfPV_, pt2PV_, rePfPV_;
  TVector3 svMinus_, svPlus_;
    
};

#endif

