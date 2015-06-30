#define hTTCPTreeAnalyzer_cxx
#include "hTTCPTreeAnalyzer.h"
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>

#include <set>
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <utility>

//declare helper functins
TVector3 impactParameter(const TVector3&, 
			 const TVector3&, 
			 const TLorentzVector&);
std::pair<float,float> angleBetweenPlanes(const TLorentzVector&, const TLorentzVector&,
					  const TLorentzVector&, const TLorentzVector&);
float deltaPhi(float, float);
float deltaR2(float, float, float, float);
float deltaR2(const TLorentzVector&, const TLorentzVector&);
enum tauDecayModes {kElectron, kMuon, 
		    kOneProng0pi0, kOneProng1pi0, kOneProng2pi0, kOneProng3pi0,
		    kThreeProng0pi0, kThreeProng1pi0,
		    kOther, kUndefined};

bool isOneProng(int decMode){
  if(decMode==kOneProng0pi0 ||
     decMode==kOneProng1pi0 ||
     decMode==kOneProng2pi0 ||
     decMode==kOneProng3pi0 ) return true;
  else return false;
}

bool isLepton(int decMode){
  if(decMode==kElectron || decMode==kMuon) return true;
  else return false;
}

Int_t hTTCPTreeAnalyzer::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  ////
  // Minimal kinematic selection
  // lepton+lepton
  if( isLepton(decModeMinus) && !isLepton(decModePlus) ){
    if(visTauMinus->Pt()>20 &&
       std::abs(visTauMinus->Eta())<2.4 &&
       visTauPlus->Pt()>20 &&
       std::abs(visTauPlus->Eta())<2.4 ) return 1;
    else return -1;
  }
  // lepton+tau
  if( isLepton(decModeMinus) && !isLepton(decModePlus) ){
    if(visTauMinus->Pt()>20 &&
       std::abs(visTauMinus->Eta())<2.1 &&
       visTauPlus->Pt()>25 &&
       std::abs(visTauPlus->Eta())<2.3 ) return 1;
    else return -1;
  }
  if( !isLepton(decModeMinus) && isLepton(decModePlus) ){
    if(visTauMinus->Pt()>25 &&
       std::abs(visTauMinus->Eta())<2.3 &&
       visTauPlus->Pt()>20 &&
       std::abs(visTauPlus->Eta())<2.1 ) return 1;
    else return -1;
  }
  // tau+tau
  if( !isLepton(decModeMinus) && !isLepton(decModePlus) ){
    if(visTauMinus->Pt()>40 &&
       std::abs(visTauMinus->Eta())<2.3 &&
       visTauPlus->Pt()>40 &&
       std::abs(visTauPlus->Eta())<2.3 ) return 1;
    else return -1;
  }
  //undefined
  return -1;
}

void hTTCPTreeAnalyzer::Loop(Long64_t maxN, bool select)
{
//   In a ROOT session, you can do:
//      Root > .L hTTCPTreeAnalyzer.C
//      Root > hTTCPTreeAnalyzer t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   //define histograms (move to constructor?)
   TH1F *hPhi[4][3];
   hPhi[0][0] = new TH1F("hPhi0", "Angle between tau decay planes (#pi-#pi); #phi^{*} [rad]; Events",32,0,TMath::Pi());
   hPhi[1][0] = new TH1F("hPhi1", "Angle between tau decay planes (1p-1p); #phi^{*} [rad]; Events",32,0,TMath::Pi());
   hPhi[2][0] = new TH1F("hPhi2", "Angle between tau decay planes (l-#pi); #phi^{*} [rad]; Events",32,0,TMath::Pi());
   hPhi[3][0] = new TH1F("hPhi3", "Angle between tau decay planes (l-1p); #phi^{*} [rad]; Events",32,0,TMath::Pi());
   TH1F *hRho[4][3];
   hRho[0][0] = new TH1F("hRho0", "Angle between visible taus in Higgs frame (#pi-#pi); #rho^{*} [rad]; Events",32,0.5*TMath::Pi(),TMath::Pi());
   hRho[1][0] = new TH1F("hRho1", "Angle between visible taus in Higgs frame (1p-1p); #rho^{*} [rad]; Events",32,0.5*TMath::Pi(),TMath::Pi());
   hRho[2][0] = new TH1F("hRho2", "Angle between visible taus in Higgs frame (l-#pi); #rho^{*} [rad]; Events",32,0.5*TMath::Pi(),TMath::Pi());
   hRho[3][0] = new TH1F("hRho3", "Angle between visible taus in Higgs frame (l-1p); #rho^{*} [rad]; Events",32,0.5*TMath::Pi(),TMath::Pi());
   TH1F *hPhi2[4][3];
   hPhi2[0][0] = new TH1F("hPhi20", "Angle between 1p-IP planes (#pi-#pi); #phi^{*} [rad]; Events",32,0,TMath::Pi());
   hPhi2[1][0] = new TH1F("hPhi21", "Angle between 1p-IP planes (1p-1p); #phi^{*} [rad]; Events",32,0,TMath::Pi());
   hPhi2[2][0] = new TH1F("hPhi22", "Angle between 1p-IP planes (l-#pi); #phi^{*} [rad]; Events",32,0,TMath::Pi());
   hPhi2[3][0] = new TH1F("hPhi23", "Angle between 1p-IP planes (l-1p); #phi^{*} [rad]; Events",32,0,TMath::Pi());
   TH1F *hRho2[4][3]; 
   hRho2[0][0] = new TH1F("hRho20", "Angle between impact parameters in #tau_{vis}-#tau_{vis} frame (#pi-#pi); #rho^{*} [rad]; Events",32,0.5*TMath::Pi(),TMath::Pi());
   hRho2[1][0] = new TH1F("hRho21", "Angle between impact parameters in #tau_{vis}-#tau_{vis} frame (1p-1p); #rho^{*} [rad]; Events",32,0.5*TMath::Pi(),TMath::Pi());
   hRho2[2][0] = new TH1F("hRho22", "Angle between impact parameters in #tau_{vis}-#tau_{vis} frame (l-#pi); #rho^{*} [rad]; Events",32,0.5*TMath::Pi(),TMath::Pi());
   hRho2[3][0] = new TH1F("hRho23", "Angle between impact parameters in #tau_{vis}-#tau_{vis} frame (l-1p); #rho^{*} [rad]; Events",32,0.5*TMath::Pi(),TMath::Pi());
   for(int ii=1;ii<3;++ii){
     for(int jj=0;jj<4;++jj){
       hPhi[jj][ii]=(TH1F*)hPhi[jj][0]->Clone(Form("hPhi%i%i",jj,ii));
       hRho[jj][ii]=(TH1F*)hRho[jj][0]->Clone(Form("hRho%i%i",jj,ii));
       hPhi2[jj][ii]=(TH1F*)hPhi2[jj][0]->Clone(Form("hPhi2%i%i",jj,ii));
       hRho2[jj][ii]=(TH1F*)hRho2[jj][0]->Clone(Form("hRho2%i%i",jj,ii));
     }
   }
   TH1F *hPhiRho[4][3];
   hPhiRho[0][0] = new TH1F("hPhiRho0", "Angle between leading charged and neutral in #tau_{vis}-#tau_{vis} frame, y^{-}*y^{+}>0 (#rho-#rho); #phi_{#rho}^{*} [rad]; Events",32,0,TMath::Pi());
   hPhiRho[1][0] = new TH1F("hPhiRho1", "Angle between leading charged and neutral in #tau_{vis}-#tau_{vis} frame, y^{-}*y^{+}<0 (#rho-#rho); #phi_{#rho}^{*} [rad]; Events",32,0,TMath::Pi());
   hPhiRho[2][0] = new TH1F("hPhiRho2", "Angle between leading charged and neutral in #tau_{vis}-#tau_{vis} frame, y^{-}*y^{+}>0 (1p+EM-1p+EM); #phi_{#rho}^{*} [rad]; Events",32,0,TMath::Pi());
   hPhiRho[3][0] = new TH1F("hPhiRho3", "Angle between leading charged and neutral in #tau_{vis}-#tau_{vis} frame, y^{-}*y^{+}<0 (1p+EM-1p+EM); #phi_{#rho}^{*} [rad]; Events",32,0,TMath::Pi());
   TH1F *hPhiRho2[4][3];
   hPhiRho2[0][0] = new TH1F("hPhiRho20", "Angle between leading charged and neutral in #tau_{vis}-#tau_{vis} frame, y_{LAB}^{-}*y_{LAB}^{+}>0 (#rho-#rho); #phi_{#rho}^{*} [rad]; Events",32,0,TMath::Pi());
   hPhiRho2[1][0] = new TH1F("hPhiRho21", "Angle between leading charged and neutral in #tau_{vis}-#tau_{vis} frame, y_{LAB}^{-}*y_{LAB}^{+}<0 (#rho-#rho); #phi_{#rho}^{*} [rad]; Events",32,0,TMath::Pi());
   hPhiRho2[2][0] = new TH1F("hPhiRho22", "Angle between leading charged and neutral in #tau_{vis}-#tau_{vis} frame, y_{LAB}^{-}*y_{LAB}^{+}>0 (1p+EM-1p+EM); #phi_{#rho}^{*} [rad]; Events",32,0,TMath::Pi());
   hPhiRho2[3][0] = new TH1F("hPhiRho23", "Angle between leading charged and neutral in #tau_{vis}-#tau_{vis} frame, y_{LAB}^{-}*y_{LAB}^{+}<0 (1p+EM-1p+EM); #phi_{#rho}^{*} [rad]; Events",32,0,TMath::Pi());
   for(int ii=1;ii<3;++ii){
     for(int jj=0;jj<4;++jj){
       hPhiRho[jj][ii]=(TH1F*)hPhiRho[jj][0]->Clone(Form("hPhiRho%i%i",jj,ii));
       hPhiRho2[jj][ii]=(TH1F*)hPhiRho2[jj][0]->Clone(Form("hPhiRho2%i%i",jj,ii));
     }
   }
   //to test consistency
   TH1F *hDPhi = new TH1F("hDPhi", "; #Delta#phi^{*} [rad]; Events",32,-0.5*TMath::Pi(),0.5*TMath::Pi());
   TH1F *hDPhi2 = new TH1F("hDPhi2", "; #Delta#phi^{*} [rad]; Events",32,-0.5*TMath::Pi(),0.5*TMath::Pi());
   TH1F *hDPhiRho = new TH1F("hDPhiRho", "; #Delta#phi^{*} [rad]; Events",32,-0.5*TMath::Pi(),0.5*TMath::Pi());

   //Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nentries = fChain->GetEntries();
   maxN = (maxN < 0 || maxN > nentries) ? nentries : maxN;

   std::cout<<maxN<<" entries will be processed"<<std::endl;

   Long64_t nbytes = 0, nb = 0;
   Long64_t nProcessed = 0, nSelected = 0;
   for (Long64_t jentry=0; jentry<maxN;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      if(jentry%1000==0) std::cout<<"Processing entry "<<jentry<<std::endl;
      nProcessed++;
      nb = fChain->GetEntry(jentry);   nbytes += nb;      
      // if (Cut(ientry) < 0) continue;
      if( select && Cut(ientry)<0 ) continue;
      //skip 3-prong and unknown decays
      if( !( isOneProng(decModeMinus)||isLepton(decModeMinus) ) ||
	  !(isOneProng(decModePlus)||isLepton(decModePlus) ) ) continue;
      nSelected++;

      std::pair<float,float> angles, angles2, anglesRho;
      angles = angleBetweenPlanes(*tauMinus,*piMinus,*tauPlus,*piPlus);
      //      angles2 = angleBetweenPlanes(*piMinus,TLorentzVector(*nPiMinus,0),
      //				   *piPlus,TLorentzVector(*nPiPlus,0));
      TVector3 nMinus = impactParameter(*thePV,*svMinus,*piMinus);
      TVector3 nPlus  = impactParameter(*thePV,*svPlus,*piPlus);
      angles2 = angleBetweenPlanes(*piMinus,TLorentzVector(nMinus,0),
      				   *piPlus,TLorentzVector(nPlus,0));
      //angles for visible tau - leading charge particle system in visible rest frame
      anglesRho = angleBetweenPlanes(*visTauMinus, *piMinus,
				     *visTauPlus, *piPlus);
      unsigned int iBos = 0;
      switch( std::abs(bosonId) ) {
        case 25:
	  iBos = 0;
	  break;
        case 35:
	  iBos = 0;
	  break;
        case 36:
	  iBos = 1;
	  break;
        case 23:
	  iBos = 2;
	  break;
        default:
	  iBos = 0;	  
	  break;
      }
      //std::cout<<"boson: "<<bosonId<<", "<<iBos<<std::endl;

      if(decModeMinus==kOneProng0pi0 && decModePlus==kOneProng0pi0){
	hPhi[0][iBos]->Fill(angles.first); hRho[0][iBos]->Fill(angles.second);      
	hPhi2[0][iBos]->Fill(angles2.first); hRho2[0][iBos]->Fill(angles2.second);
      }
      if( isOneProng(decModeMinus) && isOneProng(decModePlus) ){
	hPhi[1][iBos]->Fill(angles.first); hRho[1][iBos]->Fill(angles.second);      
	hPhi2[1][iBos]->Fill(angles2.first); hRho2[1][iBos]->Fill(angles2.second);	
      }
      if( ( decModeMinus==kOneProng0pi0 && isLepton(decModePlus) ) ||
	  ( isLepton(decModeMinus) && decModePlus==kOneProng0pi0) ){
	hPhi[2][iBos]->Fill(angles.first); hRho[2][iBos]->Fill(angles.second);      
	hPhi2[2][iBos]->Fill(angles2.first); hRho2[2][iBos]->Fill(angles2.second);	
      }
      if( ( isOneProng(decModeMinus) && isLepton(decModePlus) ) ||
	  ( isLepton(decModeMinus) && isOneProng(decModePlus) ) ){
	hPhi[3][iBos]->Fill(angles.first); hRho[3][iBos]->Fill(angles.second);      
	hPhi2[3][iBos]->Fill(angles2.first); hRho2[3][iBos]->Fill(angles2.second);	
      }
      if(decModeMinus==kOneProng1pi0 && decModePlus==kOneProng1pi0 ){
	if(yPlus*yMinus>0)
	  hPhiRho[0][iBos]->Fill(anglesRho.first);
	else
	  hPhiRho[1][iBos]->Fill(anglesRho.first);
	if(yPlusLab2*yMinusLab2>0)
	  hPhiRho2[0][iBos]->Fill(anglesRho.first);
	else
	  hPhiRho2[1][iBos]->Fill(anglesRho.first);
      }
      if( isOneProng(decModeMinus) && decModeMinus!=kOneProng0pi0 && 
	  isOneProng(decModePlus) && decModePlus!=kOneProng0pi0 ){
	if(yPlus*yMinus>0)
	  hPhiRho[2][iBos]->Fill(anglesRho.first);
	else
	  hPhiRho[3][iBos]->Fill(anglesRho.first);
	if(yPlusLab2*yMinusLab2>0)
	  hPhiRho2[2][iBos]->Fill(anglesRho.first);
	else
	  hPhiRho2[3][iBos]->Fill(anglesRho.first);
      }
      //consistency check
      hDPhi->Fill(deltaPhi(angles.first,phi));
      hDPhi2->Fill(deltaPhi(angles2.first,phi2));
      hDPhiRho->Fill(deltaPhi(anglesRho.first,phiRho));

   }//end entry loop
   if(nProcessed>0)
     std::cout<<"Selected/Processed: "<<nSelected<<" / "<<nProcessed
	      <<" = "<<(float)nSelected/(float)nProcessed<<std::endl;
   
   //move to dedicated method
   for(int j=0; j<3; ++j){
     int iCol=kBlue;
     switch( j ) {
        case 0:
	  iCol = kBlue;
	  break;
        case 1:
	  iCol = kGreen+2;
	  break;
        case 2:
	  iCol = kRed;
	  break;
        default:
	  iCol = kBlue;	  
	  break;
      }
     for(int i=0; i<4; ++i){
       hPhi[i][j]->SetStats(false);
       hPhi[i][j]->SetLineColor(iCol);
       if(std::abs(hPhi[i][j]->Integral(0,hPhi[i][j]->GetNbinsX()+1))>0)
	  hPhi[i][j]->Scale(1./hPhi[i][j]->Integral(0,hPhi[i][j]->GetNbinsX()+1));

       hRho[i][j]->SetStats(false);
       hRho[i][j]->SetLineColor(iCol);
       if(std::abs(hRho[i][j]->Integral(0,hRho[i][j]->GetNbinsX()+1))>0)
	 hRho[i][j]->Scale(1./hRho[i][j]->Integral(0,hRho[i][j]->GetNbinsX()+1));

       hPhi2[i][j]->SetStats(false);
       hPhi2[i][j]->SetLineColor(iCol);
       if(std::abs(hPhi2[i][j]->Integral(0,hPhi2[i][j]->GetNbinsX()+1))>0)
	 hPhi2[i][j]->Scale(1./hPhi2[i][j]->Integral(0,hPhi2[i][j]->GetNbinsX()+1));

       hRho2[i][j]->SetStats(false);
       hRho2[i][j]->SetLineColor(iCol);
       if(std::abs(hRho2[i][j]->Integral(0,hRho2[i][j]->GetNbinsX()+1))>0)
	 hRho2[i][j]->Scale(1./hRho2[i][j]->Integral(0,hRho2[i][j]->GetNbinsX()+1));
     }
     for(int i=0; i<4; ++i){
       hPhiRho[i][j]->SetStats(false);
       hPhiRho[i][j]->SetLineColor(iCol);
       if(std::abs(hPhiRho[i][j]->Integral(0,hPhiRho[i][j]->GetNbinsX()+1))>0)
	 hPhiRho[i][j]->Scale(1./hPhiRho[i][j]->Integral(0,hPhiRho[i][j]->GetNbinsX()+1));

       hPhiRho2[i][j]->SetStats(false);
       hPhiRho2[i][j]->SetLineColor(iCol);
       if(std::abs(hPhiRho2[i][j]->Integral(0,hPhiRho2[i][j]->GetNbinsX()+1))>0)
	 hPhiRho2[i][j]->Scale(1./hPhiRho2[i][j]->Integral(0,hPhiRho2[i][j]->GetNbinsX()+1));

     }
   }
   std::string sel;
   if(select) sel="_sel";

   TCanvas* c1 = new TCanvas("c1","TTree test example",800,800);
   c1->Clear(); c1->Divide(1,1);
   c1->SetLogy(0);
   for(int i=0; i<4; ++i){
     hPhi[i][0]->Draw();
     hPhi[i][1]->Draw("same");
     hPhi[i][2]->Draw("same");
     c1->Print((std::string(Form("phi_%i",i))+sel+std::string(".png")).c_str());
   }
   c1->Clear(); c1->Divide(1,1);
   c1->SetLogy();
   for(int i=0; i<4; ++i){
     hRho[i][0]->Draw();
     hRho[i][1]->Draw("same");
     hRho[i][2]->Draw("same");
     c1->Print((std::string(Form("rho_%i",i))+sel+std::string(".png")).c_str());
   }
   c1->SetLogy(0);
   c1->Clear(); c1->Divide(1,1);
   for(int i=0; i<4; ++i){
     hPhi2[i][0]->Draw();
     hPhi2[i][1]->Draw("same");
     hPhi2[i][2]->Draw("same");
     c1->Print((std::string(Form("phi2_%i",i))+sel+std::string(".png")).c_str());
   }
   c1->Clear(); c1->Divide(1,1);
   c1->SetLogy();
   for(int i=0; i<4; ++i){
     hRho2[i][0]->Draw();
     hRho2[i][1]->Draw("same");
     hRho2[i][2]->Draw("same");
     c1->Print((std::string(Form("rho2_%i",i))+sel+std::string(".png")).c_str());
   }
   c1->Clear(); c1->Divide(1,2);
   c1->cd(1)->SetLogy(0);c1->cd(2)->SetLogy(0);
   c1->cd(1);
   hPhiRho[0][0]->Draw();
   hPhiRho[0][1]->Draw("same");
   hPhiRho[0][2]->Draw("same");
   c1->cd(2);
   hPhiRho[1][0]->Draw();
   hPhiRho[1][1]->Draw("same");
   hPhiRho[1][2]->Draw("same");
   c1->Print((std::string(Form("phiRho_%i",0))+sel+std::string(".png")).c_str());
   c1->Clear(); c1->Divide(1,2);
   c1->cd(1)->SetLogy(0);c1->cd(2)->SetLogy(0);
   c1->cd(1);
   hPhiRho[2][0]->Draw();
   hPhiRho[2][1]->Draw("same");
   hPhiRho[2][2]->Draw("same");
   c1->cd(2);
   hPhiRho[3][0]->Draw();
   hPhiRho[3][1]->Draw("same");
   hPhiRho[3][2]->Draw("same");
   c1->Print((std::string(Form("phiRho_%i",1))+sel+std::string(".png")).c_str());
   c1->Clear(); c1->Divide(1,2);
   c1->cd(1)->SetLogy(0);c1->cd(2)->SetLogy(0);
   c1->cd(1);
   hPhiRho2[0][0]->Draw();
   hPhiRho2[0][1]->Draw("same");
   hPhiRho2[0][2]->Draw("same");
   c1->cd(2);
   hPhiRho2[1][0]->Draw();
   hPhiRho2[1][1]->Draw("same");
   hPhiRho2[1][2]->Draw("same");
   c1->Print((std::string(Form("phiRho2_%i",0))+sel+std::string(".png")).c_str());
   c1->Clear(); c1->Divide(1,2);
   c1->cd(1)->SetLogy(0);c1->cd(2)->SetLogy(0);
   c1->cd(1);
   hPhiRho2[2][0]->Draw();
   hPhiRho2[2][1]->Draw("same");
   hPhiRho2[2][2]->Draw("same");
   c1->cd(2);
   hPhiRho2[3][0]->Draw();
   hPhiRho2[3][1]->Draw("same");
   hPhiRho2[3][2]->Draw("same");
   c1->Print((std::string(Form("phiRho2_%i",1))+sel+std::string(".png")).c_str());
   
   c1->Clear(); c1->Divide(1,1);
   c1->SetLogy(0);
   hDPhi->Draw();
   c1->Print((std::string("dPhi")+sel+std::string(".png")).c_str());
   c1->Clear(); c1->Divide(1,1);
   c1->SetLogy(0);
   hDPhi2->Draw();
   c1->Print((std::string("dPhi2")+sel+std::string(".png")).c_str());
   c1->Clear(); c1->Divide(1,1);
   c1->SetLogy(0);
   hDPhiRho->Draw();
   c1->Print((std::string("dPhiRho")+sel+std::string(".png")).c_str());

   delete c1;
   
   //delete (move to destructor)
   for(int i=0; i<4; ++i){
     for(int j=0; j<3; ++j){
       delete hPhi[i][j];
       delete hRho[i][j];
       delete hPhi2[i][j];
       delete hRho2[i][j];
     }
   }
   for(int i=0; i<4; ++i){
     for(int j=0; j<3; ++j){
       delete hPhiRho[i][j];
       delete hPhiRho2[i][j];
     }
   }
   delete hDPhi;
   delete hDPhi2;

   ////
   return;
}

/// helper functions

TVector3 impactParameter(const TVector3& pv, 
			 const TVector3& sv, 
			 const TLorentzVector& p4){
  
  TVector3 dir = (p4.Vect()).Unit();
  TVector3 n = (sv-pv) - ((sv-pv)*dir)*dir;

  return n;
}

std::pair<float,float> angleBetweenPlanes(const TLorentzVector &prod1, 
					  const TLorentzVector &prod12,
					  const TLorentzVector &prod2, 
					  const TLorentzVector &prod22){
  //Boost to correct frame
  TVector3 boost = (prod1+prod2).BoostVector();
  TLorentzVector prod1Star = prod1; prod1Star.Boost(-boost);
  TLorentzVector prod12Star = prod12; prod12Star.Boost(-boost);
  TLorentzVector prod2Star = prod2; prod2Star.Boost(-boost);
  TLorentzVector prod22Star = prod22; prod22Star.Boost(-boost);

  //define common direction and normal vectors to decay planes
  TVector3 direction = prod1Star.Vect().Unit();
  TVector3 n1 = ( direction.Cross( prod12Star.Vect() ) ).Unit(); 
  TVector3 n2 = ( direction.Cross( prod22Star.Vect() ) ).Unit(); 

  float phi=TMath::ACos(n1*n2);
  float rho=TMath::ACos( (prod12Star.Vect().Unit() )*(prod22Star.Vect().Unit() ) );

  return std::make_pair(phi,rho);
}

float deltaPhi(float phi1, float phi2) { 
  float result = phi1 - phi2;
  while( result > TMath::Pi() ) result -= 2.*TMath::Pi();
  while( result <= -TMath::Pi() ) result += 2.*TMath::Pi();
  return result;
}
float deltaR2(float phi1, float eta1, float phi2, float eta2){
  return (eta1-eta2)*(eta1-eta2)+deltaPhi(phi1,phi2)*deltaPhi(phi1,phi2);
}

float deltaR2(const TLorentzVector &p41, const TLorentzVector &p42){
  return deltaR2(p41.Phi(),p41.Eta(),p42.Phi(),p42.Eta());
 }

