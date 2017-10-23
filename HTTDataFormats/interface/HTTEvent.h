#ifndef WarsawAnalysis_HTTDataFormats_HTTEvent_h
#define WarsawAnalysis_HTTDataFormats_HTTEvent_h

#include "TLorentzVector.h"
#include "TVector3.h"

//
// class declaration
//

class DiTauData{

 public:

  DiTauData();
  ~DiTauData();

  ///Data common for generator and reconstruction levels.
  int decModeMinus_, decModePlus_;
  float dzMinus_, dzPlus_;
  TVector3 thePV_;
  TVector3 svMinus_, svPlus_;
  TVector3 nPiPlus_, nPiMinus_;

  TVector3 nPiPlusAODvx_,  nPiPlusGenvx_,  nPiPlusRefitvx_;
  TVector3 nPiMinusAODvx_,  nPiMinusGenvx_,  nPiMinusRefitvx_; 
  
  TLorentzVector p4Sum_;
  TLorentzVector piMinus_, piPlus_;
  TLorentzVector tauMinus_, tauPlus_;
  TLorentzVector visTauMinus_, visTauPlus_;

  ///Reconstruction level only data
  TVector3 pfPV_, pt2PV_, refitPfPV_, refitPfPVNoBS_;
  bool isRefit_;
  int nTracksInRefit_;
  int pfPVIndex_, pt2PVindex_;

  int nPV_;
  
  std::vector<float> pfScore_, pt2Score_;
  std::vector<float> dzVtx_;

  int isoMVAWpMinus_, isoMVAWpPlus_;

  float pt2Sum_, ptThrust_, ptBalance_;

  ///Generator level only data


  ///Reset data members to default values.
  void clear();
};


class HTTEvent {

 public:

  HTTEvent();

  ~HTTEvent();

  ///Reset data members to default values.
  void clear();
  
  float run_, lumi_, event_;
  int bosonId_;

  DiTauData genEvent_, recoEvent_;
    
};

#endif

