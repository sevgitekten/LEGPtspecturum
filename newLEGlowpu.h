//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun May  3 22:59:03 2020 by ROOT version 6.08/02
// from TTree ProcessedTree/ProcessedTree
// found on file: H_LOEG.root
//////////////////////////////////////////////////////////

#ifndef newLEGlowpu_h
#define newLEGlowpu_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "Math/GenVector/PxPyPzE4D.h"
using namespace std;
class newLEGlowpu {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
    static constexpr Int_t kMaxPFMet__et = 1;
    static constexpr Int_t kMaxPFMet__sumEt = 1;
    static constexpr Int_t kMaxPFMet__phi = 1;
    static constexpr Int_t kMaxPFMetT0__et = 1;
    static constexpr Int_t kMaxPFMetT0__sumEt = 1;
    static constexpr Int_t kMaxPFMetT0__phi = 1;
    static constexpr Int_t kMaxPFMetT0T1__et = 1;
    static constexpr Int_t kMaxPFMetT0T1__sumEt = 1;
    static constexpr Int_t kMaxPFMetT0T1__phi = 1;
    static constexpr Int_t kMaxTriggerDecision = 1;
    static constexpr Int_t kMaxFilterDecision = 1;
    static constexpr Int_t kMaxL1Prescale = 1;
    static constexpr Int_t kMaxHLTPrescale = 1;
    static constexpr Int_t kMaxHLTObj = 1;
    static constexpr Int_t kMaxL1Obj = 1;
    static constexpr Int_t kMaxGenJets_ = 1;
    static constexpr Int_t kMaxPFJetsCHS_ = 51;
    static constexpr Int_t kMaxPFJetsCHS__genIdx = 51;
    static constexpr Int_t kMaxPFJetsCHS__genR = 51;
    static constexpr Int_t kMaxPFJetsCHS__cor = 51;
    static constexpr Int_t kMaxPFJetsCHS__jecLabels = 51;
    static constexpr Int_t kMaxPFJetsCHS__unc = 51;
    static constexpr Int_t kMaxPFJetsCHS__uncSrc = 51;
    static constexpr Int_t kMaxPFJetsCHS__area = 51;
    static constexpr Int_t kMaxPFJetsCHS__looseID = 51;
    static constexpr Int_t kMaxPFJetsCHS__tightID = 51;
    static constexpr Int_t kMaxPFJetsCHS__pfCombinedCvsL = 51;
    static constexpr Int_t kMaxPFJetsCHS__pfCombinedCvsB = 51;
    static constexpr Int_t kMaxPFJetsCHS__pfDeepCSVb = 51;
    static constexpr Int_t kMaxPFJetsCHS__pfDeepCSVc = 51;
    static constexpr Int_t kMaxPFJetsCHS__pfDeepCSVl = 51;
    static constexpr Int_t kMaxPFJetsCHS__pfDeepCSVbb = 51;
    static constexpr Int_t kMaxPFJetsCHS__pfBTag_JetProb = 51;
    static constexpr Int_t kMaxPFJetsCHS__pfBTag_CombInclSecVtxV2 = 51;
    static constexpr Int_t kMaxPFJetsCHS__pfBTag_CombMVAV2 = 51;
    static constexpr Int_t kMaxPFJetsCHS__QGL = 51;
    static constexpr Int_t kMaxPFJetsCHS__QGAx2 = 51;
    static constexpr Int_t kMaxPFJetsCHS__QGMul = 51;
    static constexpr Int_t kMaxPFJetsCHS__QGPtD = 51;
    static constexpr Int_t kMaxPFJetsCHS__partonFlavour = 51;
    static constexpr Int_t kMaxPFJetsCHS__partonFlavourPhysicsDef = 51;
    static constexpr Int_t kMaxPFJetsCHS__hadronFlavour = 51;
    static constexpr Int_t kMaxPFJetsCHS__chf = 51;
    static constexpr Int_t kMaxPFJetsCHS__nhf = 51;
    static constexpr Int_t kMaxPFJetsCHS__nemf = 51;
    static constexpr Int_t kMaxPFJetsCHS__cemf = 51;
    static constexpr Int_t kMaxPFJetsCHS__muf = 51;
    static constexpr Int_t kMaxPFJetsCHS__hf_hf = 51;
    static constexpr Int_t kMaxPFJetsCHS__hf_phf = 51;
    static constexpr Int_t kMaxPFJetsCHS__hf_hm = 51;
    static constexpr Int_t kMaxPFJetsCHS__hf_phm = 51;
    static constexpr Int_t kMaxPFJetsCHS__chm = 51;
    static constexpr Int_t kMaxPFJetsCHS__nhm = 51;
    static constexpr Int_t kMaxPFJetsCHS__phm = 51;
    static constexpr Int_t kMaxPFJetsCHS__elm = 51;
    static constexpr Int_t kMaxPFJetsCHS__mum = 51;
    static constexpr Int_t kMaxPFJetsCHS__ncand = 51;
    static constexpr Int_t kMaxPFJetsCHS__cm = 51;
    static constexpr Int_t kMaxPFJetsCHS__betaPrime = 51;
    static constexpr Int_t kMaxPFJetsCHS__mpuTrk = 51;
    static constexpr Int_t kMaxPFJetsCHS__mlvTrk = 51;
    static constexpr Int_t kMaxPFJetsCHS__mjtTrk = 51;
    static constexpr Int_t kMaxPFJetsCHS__hof = 51;
    static constexpr Int_t kMaxPFJetsCHS__pujid = 51;
    static constexpr Int_t kMaxgenFlavour = 1;
    static constexpr Int_t kMaxgenFlavourHadron = 1;
    static constexpr Int_t kMaxgenFlavourPartonPhysicsDef = 1;
    static constexpr Int_t kMaxgenBPt = 1;

   // Declaration of leaf types
 //QCDEvent        *events;
   Bool_t          EvtHdr__mIsPVgood;
   Bool_t          EvtHdr__mHCALNoiseNoMinZ;
   Int_t           EvtHdr__mRun;
   Long64_t        EvtHdr__mEvent;
   Int_t           EvtHdr__mLumi;
   Int_t           EvtHdr__mBunch;
   Int_t           EvtHdr__mNVtx;
   Int_t           EvtHdr__mNVtxGood;
   Int_t           EvtHdr__mOOTPUEarly;
   Int_t           EvtHdr__mOOTPULate;
   Int_t           EvtHdr__mINTPU;
   Int_t           EvtHdr__mNBX;
   Float_t         EvtHdr__mPVndof;
   Float_t         EvtHdr__mTrPu;
   Float_t         EvtHdr__mPVx;
   Float_t         EvtHdr__mPVy;
   Float_t         EvtHdr__mPVz;
   Float_t         EvtHdr__mBSx;
   Float_t         EvtHdr__mBSy;
   Float_t         EvtHdr__mBSz;
   Float_t         EvtHdr__mPthat;
   Float_t         EvtHdr__mWeight;
   Float_t         EvtHdr__mCaloRho;
   Float_t         EvtHdr__mPFRho;
   Float_t         PFMet__et_;
   Float_t         PFMet__sumEt_;
   Float_t         PFMet__phi_;
   Float_t         PFMetT0__et_;
   Float_t         PFMetT0__sumEt_;
   Float_t         PFMetT0__phi_;
   Float_t         PFMetT0T1__et_;
   Float_t         PFMetT0T1__sumEt_;
   Float_t         PFMetT0T1__phi_;
   vector<int>     TriggerDecision_;
   vector<int>     FilterDecision_;
   vector<int>     L1Prescale_;
   vector<int>     HLTPrescale_;
 //vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > HLTObj_;
 //vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > L1Obj_;
   Int_t           GenJets__;
   Double_t        GenJets__fCoordinates_fX[kMaxGenJets_];   //[GenJets__]
   Double_t        GenJets__fCoordinates_fY[kMaxGenJets_];   //[GenJets__]
   Double_t        GenJets__fCoordinates_fZ[kMaxGenJets_];   //[GenJets__]
   Double_t        GenJets__fCoordinates_fT[kMaxGenJets_];   //[GenJets__]
   Int_t           PFJetsCHS__;
   Double_t        PFJetsCHS__P4__fCoordinates_fX[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Double_t        PFJetsCHS__P4__fCoordinates_fY[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Double_t        PFJetsCHS__P4__fCoordinates_fZ[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Double_t        PFJetsCHS__P4__fCoordinates_fT[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__genIdx_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__genR_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__cor_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   vector<double>  PFJetsCHS__jecLabels_[kMaxPFJetsCHS_];
   Float_t         PFJetsCHS__unc_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   vector<float>   PFJetsCHS__uncSrc_[kMaxPFJetsCHS_];
   Float_t         PFJetsCHS__area_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Bool_t          PFJetsCHS__looseID_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Bool_t          PFJetsCHS__tightID_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__pfCombinedCvsL_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__pfCombinedCvsB_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__pfDeepCSVb_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__pfDeepCSVc_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__pfDeepCSVl_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__pfDeepCSVbb_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__pfBTag_JetProb_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__pfBTag_CombInclSecVtxV2_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__pfBTag_CombMVAV2_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__QGL_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__QGAx2_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__QGMul_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__QGPtD_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__partonFlavour_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__partonFlavourPhysicsDef_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__hadronFlavour_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__chf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__nhf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__nemf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__cemf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__muf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__hf_hf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__hf_phf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__hf_hm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__hf_phm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__chm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__nhm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__phm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__elm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__mum_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__ncand_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__cm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__betaPrime_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__mpuTrk_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__mlvTrk_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__mjtTrk_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__hof_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__pujid_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   vector<int>     genFlavour_;
   vector<int>     genFlavourHadron_;
   vector<int>     genFlavourPartonPhysicsDef_;
   vector<float>   genBPt_;

   // List of branches
   TBranch        *b_events_EvtHdr__mIsPVgood;   //!
   TBranch        *b_events_EvtHdr__mHCALNoiseNoMinZ;   //!
   TBranch        *b_events_EvtHdr__mRun;   //!
   TBranch        *b_events_EvtHdr__mEvent;   //!
   TBranch        *b_events_EvtHdr__mLumi;   //!
   TBranch        *b_events_EvtHdr__mBunch;   //!
   TBranch        *b_events_EvtHdr__mNVtx;   //!
   TBranch        *b_events_EvtHdr__mNVtxGood;   //!
   TBranch        *b_events_EvtHdr__mOOTPUEarly;   //!
   TBranch        *b_events_EvtHdr__mOOTPULate;   //!
   TBranch        *b_events_EvtHdr__mINTPU;   //!
   TBranch        *b_events_EvtHdr__mNBX;   //!
   TBranch        *b_events_EvtHdr__mPVndof;   //!
   TBranch        *b_events_EvtHdr__mTrPu;   //!
   TBranch        *b_events_EvtHdr__mPVx;   //!
   TBranch        *b_events_EvtHdr__mPVy;   //!
   TBranch        *b_events_EvtHdr__mPVz;   //!
   TBranch        *b_events_EvtHdr__mBSx;   //!
   TBranch        *b_events_EvtHdr__mBSy;   //!
   TBranch        *b_events_EvtHdr__mBSz;   //!
   TBranch        *b_events_EvtHdr__mPthat;   //!
   TBranch        *b_events_EvtHdr__mWeight;   //!
   TBranch        *b_events_EvtHdr__mCaloRho;   //!
   TBranch        *b_events_EvtHdr__mPFRho;   //!
   TBranch        *b_events_PFMet__et_;   //!
   TBranch        *b_events_PFMet__sumEt_;   //!
   TBranch        *b_events_PFMet__phi_;   //!
   TBranch        *b_events_PFMetT0__et_;   //!
   TBranch        *b_events_PFMetT0__sumEt_;   //!
   TBranch        *b_events_PFMetT0__phi_;   //!
   TBranch        *b_events_PFMetT0T1__et_;   //!
   TBranch        *b_events_PFMetT0T1__sumEt_;   //!
   TBranch        *b_events_PFMetT0T1__phi_;   //!
   TBranch        *b_events_TriggerDecision_;   //!
   TBranch        *b_events_FilterDecision_;   //!
   TBranch        *b_events_L1Prescale_;   //!
   TBranch        *b_events_HLTPrescale_;   //!
   TBranch        *b_events_GenJets__;   //!
   TBranch        *b_GenJets__fCoordinates_fX;   //!
   TBranch        *b_GenJets__fCoordinates_fY;   //!
   TBranch        *b_GenJets__fCoordinates_fZ;   //!
   TBranch        *b_GenJets__fCoordinates_fT;   //!
   TBranch        *b_events_PFJetsCHS__;   //!
   TBranch        *b_PFJetsCHS__P4__fCoordinates_fX;   //!
   TBranch        *b_PFJetsCHS__P4__fCoordinates_fY;   //!
   TBranch        *b_PFJetsCHS__P4__fCoordinates_fZ;   //!
   TBranch        *b_PFJetsCHS__P4__fCoordinates_fT;   //!
   TBranch        *b_PFJetsCHS__genIdx_;   //!
   TBranch        *b_PFJetsCHS__genR_;   //!
   TBranch        *b_PFJetsCHS__cor_;   //!
   TBranch        *b_PFJetsCHS__jecLabels_;   //!
   TBranch        *b_PFJetsCHS__unc_;   //!
   TBranch        *b_PFJetsCHS__uncSrc_;   //!
   TBranch        *b_PFJetsCHS__area_;   //!
   TBranch        *b_PFJetsCHS__looseID_;   //!
   TBranch        *b_PFJetsCHS__tightID_;   //!
   TBranch        *b_PFJetsCHS__pfCombinedCvsL_;   //!
   TBranch        *b_PFJetsCHS__pfCombinedCvsB_;   //!
   TBranch        *b_PFJetsCHS__pfDeepCSVb_;   //!
   TBranch        *b_PFJetsCHS__pfDeepCSVc_;   //!
   TBranch        *b_PFJetsCHS__pfDeepCSVl_;   //!
   TBranch        *b_PFJetsCHS__pfDeepCSVbb_;   //!
   TBranch        *b_PFJetsCHS__pfBTag_JetProb_;   //!
   TBranch        *b_PFJetsCHS__pfBTag_CombInclSecVtxV2_;   //!
   TBranch        *b_PFJetsCHS__pfBTag_CombMVAV2_;   //!
   TBranch        *b_PFJetsCHS__QGL_;   //!
   TBranch        *b_PFJetsCHS__QGAx2_;   //!
   TBranch        *b_PFJetsCHS__QGMul_;   //!
   TBranch        *b_PFJetsCHS__QGPtD_;   //!
   TBranch        *b_PFJetsCHS__partonFlavour_;   //!
   TBranch        *b_PFJetsCHS__partonFlavourPhysicsDef_;   //!
   TBranch        *b_PFJetsCHS__hadronFlavour_;   //!
   TBranch        *b_PFJetsCHS__chf_;   //!
   TBranch        *b_PFJetsCHS__nhf_;   //!
   TBranch        *b_PFJetsCHS__nemf_;   //!
   TBranch        *b_PFJetsCHS__cemf_;   //!
   TBranch        *b_PFJetsCHS__muf_;   //!
   TBranch        *b_PFJetsCHS__hf_hf_;   //!
   TBranch        *b_PFJetsCHS__hf_phf_;   //!
   TBranch        *b_PFJetsCHS__hf_hm_;   //!
   TBranch        *b_PFJetsCHS__hf_phm_;   //!
   TBranch        *b_PFJetsCHS__chm_;   //!
   TBranch        *b_PFJetsCHS__nhm_;   //!
   TBranch        *b_PFJetsCHS__phm_;   //!
   TBranch        *b_PFJetsCHS__elm_;   //!
   TBranch        *b_PFJetsCHS__mum_;   //!
   TBranch        *b_PFJetsCHS__ncand_;   //!
   TBranch        *b_PFJetsCHS__cm_;   //!
   TBranch        *b_PFJetsCHS__betaPrime_;   //!
   TBranch        *b_PFJetsCHS__mpuTrk_;   //!
   TBranch        *b_PFJetsCHS__mlvTrk_;   //!
   TBranch        *b_PFJetsCHS__mjtTrk_;   //!
   TBranch        *b_PFJetsCHS__hof_;   //!
   TBranch        *b_PFJetsCHS__pujid_;   //!
   TBranch        *b_events_genFlavour_;   //!
   TBranch        *b_events_genFlavourHadron_;   //!
   TBranch        *b_events_genFlavourPartonPhysicsDef_;   //!
   TBranch        *b_events_genBPt_;   //!

   newLEGlowpu(TTree *tree=0);
   virtual ~newLEGlowpu();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef newLEGlowpu_cxx
newLEGlowpu::newLEGlowpu(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("H_LOEG.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("H_LOEG.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("H_LOEG.root:/ak4");
      dir->GetObject("ProcessedTree",tree);

   }
   Init(tree);
}

newLEGlowpu::~newLEGlowpu()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t newLEGlowpu::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t newLEGlowpu::LoadTree(Long64_t entry)
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

void newLEGlowpu::Init(TTree *tree)
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

   fChain->SetBranchAddress("EvtHdr_.mIsPVgood", &EvtHdr__mIsPVgood, &b_events_EvtHdr__mIsPVgood);
   fChain->SetBranchAddress("EvtHdr_.mHCALNoiseNoMinZ", &EvtHdr__mHCALNoiseNoMinZ, &b_events_EvtHdr__mHCALNoiseNoMinZ);
   fChain->SetBranchAddress("EvtHdr_.mRun", &EvtHdr__mRun, &b_events_EvtHdr__mRun);
   fChain->SetBranchAddress("EvtHdr_.mEvent", &EvtHdr__mEvent, &b_events_EvtHdr__mEvent);
   fChain->SetBranchAddress("EvtHdr_.mLumi", &EvtHdr__mLumi, &b_events_EvtHdr__mLumi);
   fChain->SetBranchAddress("EvtHdr_.mBunch", &EvtHdr__mBunch, &b_events_EvtHdr__mBunch);
   fChain->SetBranchAddress("EvtHdr_.mNVtx", &EvtHdr__mNVtx, &b_events_EvtHdr__mNVtx);
   fChain->SetBranchAddress("EvtHdr_.mNVtxGood", &EvtHdr__mNVtxGood, &b_events_EvtHdr__mNVtxGood);
   fChain->SetBranchAddress("EvtHdr_.mOOTPUEarly", &EvtHdr__mOOTPUEarly, &b_events_EvtHdr__mOOTPUEarly);
   fChain->SetBranchAddress("EvtHdr_.mOOTPULate", &EvtHdr__mOOTPULate, &b_events_EvtHdr__mOOTPULate);
   fChain->SetBranchAddress("EvtHdr_.mINTPU", &EvtHdr__mINTPU, &b_events_EvtHdr__mINTPU);
   fChain->SetBranchAddress("EvtHdr_.mNBX", &EvtHdr__mNBX, &b_events_EvtHdr__mNBX);
   fChain->SetBranchAddress("EvtHdr_.mPVndof", &EvtHdr__mPVndof, &b_events_EvtHdr__mPVndof);
   fChain->SetBranchAddress("EvtHdr_.mTrPu", &EvtHdr__mTrPu, &b_events_EvtHdr__mTrPu);
   fChain->SetBranchAddress("EvtHdr_.mPVx", &EvtHdr__mPVx, &b_events_EvtHdr__mPVx);
   fChain->SetBranchAddress("EvtHdr_.mPVy", &EvtHdr__mPVy, &b_events_EvtHdr__mPVy);
   fChain->SetBranchAddress("EvtHdr_.mPVz", &EvtHdr__mPVz, &b_events_EvtHdr__mPVz);
   fChain->SetBranchAddress("EvtHdr_.mBSx", &EvtHdr__mBSx, &b_events_EvtHdr__mBSx);
   fChain->SetBranchAddress("EvtHdr_.mBSy", &EvtHdr__mBSy, &b_events_EvtHdr__mBSy);
   fChain->SetBranchAddress("EvtHdr_.mBSz", &EvtHdr__mBSz, &b_events_EvtHdr__mBSz);
   fChain->SetBranchAddress("EvtHdr_.mPthat", &EvtHdr__mPthat, &b_events_EvtHdr__mPthat);
   fChain->SetBranchAddress("EvtHdr_.mWeight", &EvtHdr__mWeight, &b_events_EvtHdr__mWeight);
   fChain->SetBranchAddress("EvtHdr_.mCaloRho", &EvtHdr__mCaloRho, &b_events_EvtHdr__mCaloRho);
   fChain->SetBranchAddress("EvtHdr_.mPFRho", &EvtHdr__mPFRho, &b_events_EvtHdr__mPFRho);
   fChain->SetBranchAddress("PFMet_.et_", &PFMet__et_, &b_events_PFMet__et_);
   fChain->SetBranchAddress("PFMet_.sumEt_", &PFMet__sumEt_, &b_events_PFMet__sumEt_);
   fChain->SetBranchAddress("PFMet_.phi_", &PFMet__phi_, &b_events_PFMet__phi_);
   fChain->SetBranchAddress("PFMetT0_.et_", &PFMetT0__et_, &b_events_PFMetT0__et_);
   fChain->SetBranchAddress("PFMetT0_.sumEt_", &PFMetT0__sumEt_, &b_events_PFMetT0__sumEt_);
   fChain->SetBranchAddress("PFMetT0_.phi_", &PFMetT0__phi_, &b_events_PFMetT0__phi_);
   fChain->SetBranchAddress("PFMetT0T1_.et_", &PFMetT0T1__et_, &b_events_PFMetT0T1__et_);
   fChain->SetBranchAddress("PFMetT0T1_.sumEt_", &PFMetT0T1__sumEt_, &b_events_PFMetT0T1__sumEt_);
   fChain->SetBranchAddress("PFMetT0T1_.phi_", &PFMetT0T1__phi_, &b_events_PFMetT0T1__phi_);
   fChain->SetBranchAddress("TriggerDecision_", &TriggerDecision_, &b_events_TriggerDecision_);
   fChain->SetBranchAddress("FilterDecision_", &FilterDecision_, &b_events_FilterDecision_);
   fChain->SetBranchAddress("L1Prescale_", &L1Prescale_, &b_events_L1Prescale_);
   fChain->SetBranchAddress("HLTPrescale_", &HLTPrescale_, &b_events_HLTPrescale_);
   fChain->SetBranchAddress("GenJets_", &GenJets__, &b_events_GenJets__);
   fChain->SetBranchAddress("GenJets_.fCoordinates.fX", &GenJets__fCoordinates_fX, &b_GenJets__fCoordinates_fX);
   fChain->SetBranchAddress("GenJets_.fCoordinates.fY", &GenJets__fCoordinates_fY, &b_GenJets__fCoordinates_fY);
   fChain->SetBranchAddress("GenJets_.fCoordinates.fZ", &GenJets__fCoordinates_fZ, &b_GenJets__fCoordinates_fZ);
   fChain->SetBranchAddress("GenJets_.fCoordinates.fT", &GenJets__fCoordinates_fT, &b_GenJets__fCoordinates_fT);
   fChain->SetBranchAddress("PFJetsCHS_", &PFJetsCHS__, &b_events_PFJetsCHS__);
   fChain->SetBranchAddress("PFJetsCHS_.P4_.fCoordinates.fX", PFJetsCHS__P4__fCoordinates_fX, &b_PFJetsCHS__P4__fCoordinates_fX);
   fChain->SetBranchAddress("PFJetsCHS_.P4_.fCoordinates.fY", PFJetsCHS__P4__fCoordinates_fY, &b_PFJetsCHS__P4__fCoordinates_fY);
   fChain->SetBranchAddress("PFJetsCHS_.P4_.fCoordinates.fZ", PFJetsCHS__P4__fCoordinates_fZ, &b_PFJetsCHS__P4__fCoordinates_fZ);
   fChain->SetBranchAddress("PFJetsCHS_.P4_.fCoordinates.fT", PFJetsCHS__P4__fCoordinates_fT, &b_PFJetsCHS__P4__fCoordinates_fT);
   fChain->SetBranchAddress("PFJetsCHS_.genIdx_", PFJetsCHS__genIdx_, &b_PFJetsCHS__genIdx_);
   fChain->SetBranchAddress("PFJetsCHS_.genR_", PFJetsCHS__genR_, &b_PFJetsCHS__genR_);
   fChain->SetBranchAddress("PFJetsCHS_.cor_", PFJetsCHS__cor_, &b_PFJetsCHS__cor_);
   fChain->SetBranchAddress("PFJetsCHS_.jecLabels_", PFJetsCHS__jecLabels_, &b_PFJetsCHS__jecLabels_);
   fChain->SetBranchAddress("PFJetsCHS_.unc_", PFJetsCHS__unc_, &b_PFJetsCHS__unc_);
   fChain->SetBranchAddress("PFJetsCHS_.uncSrc_", PFJetsCHS__uncSrc_, &b_PFJetsCHS__uncSrc_);
   fChain->SetBranchAddress("PFJetsCHS_.area_", PFJetsCHS__area_, &b_PFJetsCHS__area_);
   fChain->SetBranchAddress("PFJetsCHS_.looseID_", PFJetsCHS__looseID_, &b_PFJetsCHS__looseID_);
   fChain->SetBranchAddress("PFJetsCHS_.tightID_", PFJetsCHS__tightID_, &b_PFJetsCHS__tightID_);
   fChain->SetBranchAddress("PFJetsCHS_.pfCombinedCvsL_", PFJetsCHS__pfCombinedCvsL_, &b_PFJetsCHS__pfCombinedCvsL_);
   fChain->SetBranchAddress("PFJetsCHS_.pfCombinedCvsB_", PFJetsCHS__pfCombinedCvsB_, &b_PFJetsCHS__pfCombinedCvsB_);
   fChain->SetBranchAddress("PFJetsCHS_.pfDeepCSVb_", PFJetsCHS__pfDeepCSVb_, &b_PFJetsCHS__pfDeepCSVb_);
   fChain->SetBranchAddress("PFJetsCHS_.pfDeepCSVc_", PFJetsCHS__pfDeepCSVc_, &b_PFJetsCHS__pfDeepCSVc_);
   fChain->SetBranchAddress("PFJetsCHS_.pfDeepCSVl_", PFJetsCHS__pfDeepCSVl_, &b_PFJetsCHS__pfDeepCSVl_);
   fChain->SetBranchAddress("PFJetsCHS_.pfDeepCSVbb_", PFJetsCHS__pfDeepCSVbb_, &b_PFJetsCHS__pfDeepCSVbb_);
   fChain->SetBranchAddress("PFJetsCHS_.pfBTag_JetProb_", PFJetsCHS__pfBTag_JetProb_, &b_PFJetsCHS__pfBTag_JetProb_);
   fChain->SetBranchAddress("PFJetsCHS_.pfBTag_CombInclSecVtxV2_", PFJetsCHS__pfBTag_CombInclSecVtxV2_, &b_PFJetsCHS__pfBTag_CombInclSecVtxV2_);
   fChain->SetBranchAddress("PFJetsCHS_.pfBTag_CombMVAV2_", PFJetsCHS__pfBTag_CombMVAV2_, &b_PFJetsCHS__pfBTag_CombMVAV2_);
   fChain->SetBranchAddress("PFJetsCHS_.QGL_", PFJetsCHS__QGL_, &b_PFJetsCHS__QGL_);
   fChain->SetBranchAddress("PFJetsCHS_.QGAx2_", PFJetsCHS__QGAx2_, &b_PFJetsCHS__QGAx2_);
   fChain->SetBranchAddress("PFJetsCHS_.QGMul_", PFJetsCHS__QGMul_, &b_PFJetsCHS__QGMul_);
   fChain->SetBranchAddress("PFJetsCHS_.QGPtD_", PFJetsCHS__QGPtD_, &b_PFJetsCHS__QGPtD_);
   fChain->SetBranchAddress("PFJetsCHS_.partonFlavour_", PFJetsCHS__partonFlavour_, &b_PFJetsCHS__partonFlavour_);
   fChain->SetBranchAddress("PFJetsCHS_.partonFlavourPhysicsDef_", PFJetsCHS__partonFlavourPhysicsDef_, &b_PFJetsCHS__partonFlavourPhysicsDef_);
   fChain->SetBranchAddress("PFJetsCHS_.hadronFlavour_", PFJetsCHS__hadronFlavour_, &b_PFJetsCHS__hadronFlavour_);
   fChain->SetBranchAddress("PFJetsCHS_.chf_", PFJetsCHS__chf_, &b_PFJetsCHS__chf_);
   fChain->SetBranchAddress("PFJetsCHS_.nhf_", PFJetsCHS__nhf_, &b_PFJetsCHS__nhf_);
   fChain->SetBranchAddress("PFJetsCHS_.nemf_", PFJetsCHS__nemf_, &b_PFJetsCHS__nemf_);
   fChain->SetBranchAddress("PFJetsCHS_.cemf_", PFJetsCHS__cemf_, &b_PFJetsCHS__cemf_);
   fChain->SetBranchAddress("PFJetsCHS_.muf_", PFJetsCHS__muf_, &b_PFJetsCHS__muf_);
   fChain->SetBranchAddress("PFJetsCHS_.hf_hf_", PFJetsCHS__hf_hf_, &b_PFJetsCHS__hf_hf_);
   fChain->SetBranchAddress("PFJetsCHS_.hf_phf_", PFJetsCHS__hf_phf_, &b_PFJetsCHS__hf_phf_);
   fChain->SetBranchAddress("PFJetsCHS_.hf_hm_", PFJetsCHS__hf_hm_, &b_PFJetsCHS__hf_hm_);
   fChain->SetBranchAddress("PFJetsCHS_.hf_phm_", PFJetsCHS__hf_phm_, &b_PFJetsCHS__hf_phm_);
   fChain->SetBranchAddress("PFJetsCHS_.chm_", PFJetsCHS__chm_, &b_PFJetsCHS__chm_);
   fChain->SetBranchAddress("PFJetsCHS_.nhm_", PFJetsCHS__nhm_, &b_PFJetsCHS__nhm_);
   fChain->SetBranchAddress("PFJetsCHS_.phm_", PFJetsCHS__phm_, &b_PFJetsCHS__phm_);
   fChain->SetBranchAddress("PFJetsCHS_.elm_", PFJetsCHS__elm_, &b_PFJetsCHS__elm_);
   fChain->SetBranchAddress("PFJetsCHS_.mum_", PFJetsCHS__mum_, &b_PFJetsCHS__mum_);
   fChain->SetBranchAddress("PFJetsCHS_.ncand_", PFJetsCHS__ncand_, &b_PFJetsCHS__ncand_);
   fChain->SetBranchAddress("PFJetsCHS_.cm_", PFJetsCHS__cm_, &b_PFJetsCHS__cm_);
   fChain->SetBranchAddress("PFJetsCHS_.betaPrime_", PFJetsCHS__betaPrime_, &b_PFJetsCHS__betaPrime_);
   fChain->SetBranchAddress("PFJetsCHS_.mpuTrk_", PFJetsCHS__mpuTrk_, &b_PFJetsCHS__mpuTrk_);
   fChain->SetBranchAddress("PFJetsCHS_.mlvTrk_", PFJetsCHS__mlvTrk_, &b_PFJetsCHS__mlvTrk_);
   fChain->SetBranchAddress("PFJetsCHS_.mjtTrk_", PFJetsCHS__mjtTrk_, &b_PFJetsCHS__mjtTrk_);
   fChain->SetBranchAddress("PFJetsCHS_.hof_", PFJetsCHS__hof_, &b_PFJetsCHS__hof_);
   fChain->SetBranchAddress("PFJetsCHS_.pujid_", PFJetsCHS__pujid_, &b_PFJetsCHS__pujid_);
   fChain->SetBranchAddress("genFlavour_", &genFlavour_, &b_events_genFlavour_);
   fChain->SetBranchAddress("genFlavourHadron_", &genFlavourHadron_, &b_events_genFlavourHadron_);
   fChain->SetBranchAddress("genFlavourPartonPhysicsDef_", &genFlavourPartonPhysicsDef_, &b_events_genFlavourPartonPhysicsDef_);
   fChain->SetBranchAddress("genBPt_", &genBPt_, &b_events_genBPt_);
   Notify();
}

Bool_t newLEGlowpu::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void newLEGlowpu::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t newLEGlowpu::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef newLEGlowpu_cxx
