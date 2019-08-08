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
  TLorentzVector trkMinus_, trkPlus_;
  TLorentzVector scMinus_, scPlus_;
  TLorentzVector neutralMinus_, neutralPlus_;
  TLorentzVector tauMinus_, tauPlus_;
  TLorentzVector visTauMinus_, visTauPlus_;
  TLorentzVector pi0Minus_, pi0Plus_;

  ///Reconstruction level only data
  TVector3 pfPV_, pt2PV_, refitPfPV_, refitPfPVNoBS_;
  bool isRefit_;
  int nTracksInRefit_;
  int pfPVIndex_, pt2PVindex_;
  int leadIdMinus_, leadIdPlus_;
  int nGammaMinus_, nGammaPlus_;
  int nGammaInConeMinus_, nGammaInConePlus_;
  int nStripMinus_, nStripPlus_;

  std::vector<float> pfScore_, pt2Score_;
  std::vector<float> dzVtx_;

  float isoMinus_, isoPlus_;
  int isoMVAWpMinus_, isoMVAWpPlus_;
  int antiEWpMinus_, antiEWpPlus_;
  int antiE2WpMinus_, antiE2WpPlus_;
  int antiMuWpMinus_, antiMuWpPlus_;
  int matchedMinus_, matchedPlus_;
  float gammaPtSumInScMinus_, gammaPtSumInScPlus_;
  float gammaPtSumOutScMinus_, gammaPtSumOutScPlus_;
  float maxGammaPtMinus_, maxGammaPtPlus_;
  float dR2midMinus_, dR2midPlus_;

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
  std::vector<float> tauSpinnerWeight_;

  int nPV_;
  bool diMuonVeto_, diEleVeto_;

  DiTauData genEvent_, recoEvent_;
    
};

#endif

