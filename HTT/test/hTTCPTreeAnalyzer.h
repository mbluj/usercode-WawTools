//////////////////////////////////////////////////////////

#ifndef hTTCPTreeAnalyzer_h
#define hTTCPTreeAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <TLorentzVector.h>
#include <TObject.h>
#include <TVector3.h>

// Fixed size dimensions of array or collections stored in the TTree if any.

class hTTCPTreeAnalyzer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   TLorentzVector  *p4Sum;
   TLorentzVector  *metNu, *met;
   TLorentzVector  *piMinus, *piPlus;
   TLorentzVector  *tauMinus, *tauPlus;
   TLorentzVector  *visTauMinus, *visTauPlus;
   Int_t           bosonId, decModeMinus, decModePlus;
   TVector3        *thePV, *svMinus, *svPlus;
   TVector3        *nPiMinus, *nPiPlus;
   Float_t         phi,rho, phi2, phiRho;
   Float_t         yMinus, yPlus, yMinus2, yPlus2;
   Float_t         yMinusLab, yPlusLab, yMinusLab2, yPlusLab2;

   // List of branches
   TBranch        *b_p4Sum, *b_metNu, *b_met;
   TBranch        *b_piMinus, *b_piPlus;
   TBranch        *b_tauMinus, *b_tauPlus;
   TBranch        *b_visTauMinus, *b_visTauPlus;
   TBranch        *b_bosonId, *b_decModeMinus, *b_decModePlus;
   TBranch        *b_thePV, *b_svMinus, *b_svPlus;
   TBranch        *b_nPiMinus, *b_nPiPlus; 
   TBranch        *b_phi, *b_rho, *b_phi2, *b_phiRho;
   TBranch        *b_yMinus, *b_yPlus, *b_yMinus2, *b_yPlus2;
   TBranch        *b_yMinusLab, *b_yPlusLab, *b_yMinusLab2, *b_yPlusLab2;

   hTTCPTreeAnalyzer(TTree *tree=0);
   virtual ~hTTCPTreeAnalyzer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Long64_t maxN = -1, bool select=false);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef hTTCPTreeAnalyzer_cxx
hTTCPTreeAnalyzer::hTTCPTreeAnalyzer(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("hTTCPTree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("hTTCPTree","hTTCPTree");
      chain->Add("hTTCPNtuple_h0.root/hTTCPTree");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

hTTCPTreeAnalyzer::~hTTCPTreeAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t hTTCPTreeAnalyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t hTTCPTreeAnalyzer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void hTTCPTreeAnalyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   //fChain->SetMakeClass(1);//MB switch off split mode to retrieve full classes

   fChain->SetBranchAddress("p4Sum.", &p4Sum, &b_p4Sum);
   fChain->SetBranchAddress("metNu.", &metNu, &b_metNu);
   fChain->SetBranchAddress("met.", &met, &b_met);
   fChain->SetBranchAddress("piMinus.", &piMinus, &b_piMinus);
   fChain->SetBranchAddress("piPlus.", &piPlus, &b_piPlus);
   fChain->SetBranchAddress("tauMinus.", &tauMinus, &b_tauMinus);
   fChain->SetBranchAddress("tauPlus.", &tauPlus, &b_tauPlus);
   fChain->SetBranchAddress("visTauMinus.", &visTauMinus, &b_visTauMinus);
   fChain->SetBranchAddress("visTauPlus.", &visTauPlus, &b_visTauPlus);
   //
   fChain->SetBranchAddress("bosonId", &bosonId, &b_bosonId);
   fChain->SetBranchAddress("decModeMinus", &decModeMinus, &b_decModeMinus);
   fChain->SetBranchAddress("decModePlus", &decModePlus, &b_decModePlus);
   //
   fChain->SetBranchAddress("thePV.", &thePV, &b_thePV);
   fChain->SetBranchAddress("svMinus.", &svMinus, &b_svMinus);
   fChain->SetBranchAddress("svPlus.", &svPlus, &b_svPlus);
   fChain->SetBranchAddress("nPiMinus.", &nPiMinus, &b_nPiMinus);
   fChain->SetBranchAddress("nPiPlus.", &nPiPlus, &b_nPiPlus);
   //
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("phi2", &phi2, &b_phi2);
   fChain->SetBranchAddress("phiRho", &phiRho, &b_phiRho);
   //
   fChain->SetBranchAddress("yMinus", &yMinus, &b_yMinus);
   fChain->SetBranchAddress("yPlus", &yPlus, &b_yPlus);
   fChain->SetBranchAddress("yMinus2", &yMinus2, &b_yMinus2);
   fChain->SetBranchAddress("yPlus2", &yPlus2, &b_yPlus2);
   fChain->SetBranchAddress("yMinusLab", &yMinusLab, &b_yMinusLab);
   fChain->SetBranchAddress("yPlusLab", &yPlusLab, &b_yPlusLab);
   fChain->SetBranchAddress("yMinusLab2", &yMinusLab2, &b_yMinusLab2);
   fChain->SetBranchAddress("yPlusLab2", &yPlusLab2, &b_yPlusLab2);

   Notify();
}

Bool_t hTTCPTreeAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void hTTCPTreeAnalyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
/* implementation moved to .C 
Int_t hTTCPTreeAnalyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
*/
#endif // #ifdef hTTCPTreeAnalyzer_cxx
