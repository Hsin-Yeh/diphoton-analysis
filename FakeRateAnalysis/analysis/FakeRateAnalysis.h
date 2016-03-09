//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar  8 09:30:18 2016 by ROOT version 6.02/13
// from TTree fTree/PhotonTree
// found on file: root://cmsxrootd.fnal.gov//store/user/abuccill/DiPhotonAnalysis/FakeRateMerged/diphoton_fakeRate_JetHT_Run2015C_25ns-16Dec2015-v1_MINIAOD_merged.root
//////////////////////////////////////////////////////////

#ifndef FakeRateAnalysis_h
#define FakeRateAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class FakeRateAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         Event_pthat;
   Float_t         Event_alphaqcd;
   Float_t         Event_alphaqed;
   Float_t         Event_qscale;
   Float_t         Event_weight;
   Int_t           Event_run;
   Int_t           Event_LS;
   Int_t           Event_evnum;
   Int_t           Event_processid;
   Int_t           Event_bx;
   Int_t           Event_orbit;
   Int_t           Event_interactingParton1PdgId;
   Int_t           Event_interactingParton2PdgId;
   Float_t         Jet_jetHT;
   Int_t           Jet_nJets;
   Double_t        Photon_pt;
   Double_t        Photon_eta;
   Double_t        Photon_phi;
   Double_t        Photon_scEta;
   Double_t        Photon_scPhi;
   Double_t        Photon_rho;
   Double_t        Photon_chargedHadIso03;
   Double_t        Photon_neutralHadIso03;
   Double_t        Photon_photonIso03;
   Double_t        Photon_rhoCorChargedHadIso03;
   Double_t        Photon_rhoCorNeutralHadIso03;
   Double_t        Photon_rhoCorPhotonIso03;
   Double_t        Photon_corPhotonIso03;
   Double_t        Photon_hadTowerOverEm;
   Double_t        Photon_r9;
   Double_t        Photon_sigmaIetaIeta;
   Double_t        Photon_sigmaIetaIeta5x5;
   Double_t        Photon_sigmaEtaEta;
   Double_t        Photon_sigmaIphiIphi;
   Double_t        Photon_sigmaPhiPhi;
   Double_t        Photon_maxEnergyXtal;
   Double_t        Photon_alphaHighPtID;
   Double_t        Photon_kappaHighPtID;
   Double_t        Photon_phoEAHighPtID;
   Double_t        Photon_chEAegmID;
   Double_t        Photon_nhEAegmID;
   Double_t        Photon_phoEAegmID;
   Bool_t          Photon_passEGMLooseID;
   Bool_t          Photon_passEGMMediumID;
   Bool_t          Photon_passEGMTightID;
   Bool_t          Photon_isEB;
   Bool_t          Photon_isEE;
   Bool_t          Photon_isEBEtaGap;
   Bool_t          Photon_isEBPhiGap;
   Bool_t          Photon_isEERingGap;
   Bool_t          Photon_isEEDeeGap;
   Bool_t          Photon_isEBEEGap;
   Bool_t          Photon_passElectronVeto;
   Bool_t          Photon_passHTowOverE;
   Bool_t          Photon_passChIso;
   Bool_t          Photon_passCorPhoIso;
   Bool_t          Photon_passSieie;
   Bool_t          Photon_passHighPtID;
   Bool_t          Photon_passChIsoDenom;
   Bool_t          Photon_passCorPhoIsoDenom;
   Bool_t          Photon_isFakeable;
   Bool_t          Photon_isNumeratorObj;
   Bool_t          Photon_isDenominatorObj;
   Bool_t          Photon_isSaturated;

   // List of branches
   TBranch        *b_Event;   //!
   TBranch        *b_Jet;   //!
   TBranch        *b_Photon;   //!

   FakeRateAnalysis(TTree *tree=0);
   virtual ~FakeRateAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef FakeRateAnalysis_cxx
FakeRateAnalysis::FakeRateAnalysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://cmsxrootd.fnal.gov//store/user/abuccill/DiPhotonAnalysis/FakeRateMerged/diphoton_fakeRate_JetHT_Run2015C_25ns-16Dec2015-v1_MINIAOD_merged.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://cmsxrootd.fnal.gov//store/user/abuccill/DiPhotonAnalysis/FakeRateMerged/diphoton_fakeRate_JetHT_Run2015C_25ns-16Dec2015-v1_MINIAOD_merged.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("root://cmsxrootd.fnal.gov//store/user/abuccill/DiPhotonAnalysis/FakeRateMerged/diphoton_fakeRate_JetHT_Run2015C_25ns-16Dec2015-v1_MINIAOD_merged.root:/diphoton");
      dir->GetObject("fTree",tree);

   }
   Init(tree);
}

FakeRateAnalysis::~FakeRateAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t FakeRateAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t FakeRateAnalysis::LoadTree(Long64_t entry)
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

void FakeRateAnalysis::Init(TTree *tree)
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
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event", &Event_pthat, &b_Event);
   fChain->SetBranchAddress("Jet", &Jet_jetHT, &b_Jet);
   fChain->SetBranchAddress("Photon", &Photon_pt, &b_Photon);
   Notify();
}

Bool_t FakeRateAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void FakeRateAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t FakeRateAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef FakeRateAnalysis_cxx