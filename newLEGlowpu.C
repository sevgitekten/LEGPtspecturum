#include <iostream>
#define newLEGlowpu_cxx
#include "newLEGlowpu.h"
//#include <fstream>
//#include <sstream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TProfile.h"
#include <set>
#include "TRandom3.h"
#include <vector>
#include "basic.h"
#include <string>
#include <map>
#include <utility>
#include "TString.h"
#include "TApplication.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include <time.h>

using namespace std;

void newLEGlowpu::Loop()
{
    //ofstream runfile;
    //runfile.open ("LEGJet25run.txt");
    //time_t start,end; double dif;  time (&start);
    
    if (fChain == 0) return; //MakeClass
    Long64_t nentries = fChain->GetEntriesFast();//MakeClass
    Long64_t nbytes = 0, nb = 0; //MakeClass
    
    cout<<"nentries"<<nentries<<endl;
    //nentries = 1000000;
    
    TFile myFile("6Mayis_LEG_17Nov2017Bv6_JEC_mikko_withoutfilter.root", "RECREATE");
  
    //CondFormat'ın icinden cektigimiz text dosyaları,degistirdiklerim var!!///
    JetCorrectorParameters *pfchs_l1 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Fall17_17Nov2017B_V6_DATA_L1FastJet_AK4PFchs.txt");
    JetCorrectorParameters *pfchs_l2 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Fall17_17Nov2017B_V6_DATA_L2Relative_AK4PFchs.txt");
    JetCorrectorParameters *pfchs_l3 = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Fall17_17Nov2017B_V6_DATA_L3Absolute_AK4PFchs.txt");
    JetCorrectorParameters *pfchs_l2l3res = new JetCorrectorParameters("CondFormats/JetMETObjects/data/Fall17_17Nov2017B_V6_DATA_L2L3Residual_AK4PFchs.txt");
    
    vector<JetCorrectorParameters> vParam_pfchs;
    vParam_pfchs.push_back(*pfchs_l1);
    vParam_pfchs.push_back(*pfchs_l2);
    vParam_pfchs.push_back(*pfchs_l3);
    vParam_pfchs.push_back(*pfchs_l2l3res);
    FactorizedJetCorrector *pfchs_jec = new FactorizedJetCorrector(vParam_pfchs);
    
  
    static const int netabins = 8;
    static const double etabins[netabins+1] = {0,0.5,1.0,1.5,2.0,2.5,3.0,3.2,4.7};
    //static const double etabins[netabins+1] = {0,0.5,1.0,1.5,2.0,2.5,3.0,3.2,3.7,4.2,4.7};
    //static const double etabins[netabins+1] = {2.5,3.0,3.2,4.7};
    const double x[8][65]=
    {
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3832, 6076, 6389}, // Eta_0.0-0.5
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3637, 5220, 5492}, // Eta_0.5-1.0
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2941, 3832}, // Eta_1.0-1.5
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2500, 2640}, // Eta_1.5-2.0
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684}, // Eta_2.0-2.5
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032}, // Eta_2.5-3.0
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032}, // Eta_3.0-3.5
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032},// Eta_3.5-4.0
        /*{10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032}, // Eta_4.0-4.5
         {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032} // Eta_4.5-5.0
         */
    };
    
    //const int nx[3] = {34,34,34};
    const int nx[8] = {64,63,58,54,48,40,40,40};
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    vector<TH1F*> vectptt[4];
    
    //TH1F *hpt[netabins];
   
    for ( int k = 0; k<netabins; ++k)
    {
        /*TH1F *hlt_0   = new TH1F(Form("hlt_0_%d",k),vectrigname[0].c_str(),nx[k],&x[k][0]);
         vectptt[0].push_back (hlt_0); /// trigger histolarının oldugu vektor
         hlt_0->Sumw2();*/
        
        TH1F *hlt_1   = new TH1F(Form("hlt_1_%d",k),"cross section",nx[k],&x[k][0]);
        vectptt[1].push_back (hlt_1);
        hlt_1->Sumw2();
        
        TH1F *hlt_2   = new TH1F(Form("hlt_2_%d",k),"cross section",nx[k],&x[k][0]);
        vectptt[2].push_back (hlt_2);
        hlt_2->Sumw2();
        
        TH1F *hlt_3   = new TH1F(Form("hlt_3_%d",k),"cross section",nx[k],&x[k][0]);
        vectptt[3].push_back (hlt_3);
        hlt_3->Sumw2();
        
    }

    TLorentzVector p4pf;
    TLorentzVector p4pf2;
    vector<float>pfchs_ptcorr;
    
    
    //trigger PT Ranges
    vector<double>vecttrigcuts={49,49,64,97,133,196,300}; //10Mart
    vector<double>vecttrigcuts1={49,49,64,97,133,196,300};
    vector<double>vecttrigcuts2={56,56,74,97,133,196,300};
    vector<double>vecttrigcuts3={43,43,64,97,133,196,300};
    //================luminosity values
    vector<double>vectlumi={12410.999433624,74815.505012496,476831.591488795,2138811.439403298,7142065.860181809,211511974.809555233}; //Subat hesaplari
    
    //=================Define the variables for Correction========================
    
    float pf_jtpt[100],pf_jt_eta[100];	//corr definition (disarda kalailir)
    
    //---------------------------events dongusu icine girdik---------
    /*Long64_t dummy=270000000; int counter=0;
    TH1F *runnocount = new TH1F("runnocount","runnocount",301,306799.5,307100.5);
    TH1F *runnocount2 = new TH1F("runnocount2","runnocount2",301,306799.5,307100.5);
    */
    for (Long64_t jentry=0; jentry<nentries;jentry++)   //MakeClass
        {
        Long64_t ientry = LoadTree(jentry);             //MakeClass
        if (ientry < 0) break;                          //MakeClass
        nb = fChain->GetEntry(jentry);   nbytes += nb;  //MakeClass
        cout<<"\r"<<"event number: "<<jentry<<"/ "<<nentries<<flush;
        
        int pf_njt =PFJetsCHS__;  //==========corr definition
        pfchs_ptcorr.resize(pf_njt);	//======corr definition
        
        if (TriggerDecision_.empty()) continue;
        //if (FilterDecision_.size()!=0)  continue;
        
        //if(counter==0) dummy=EvtHdr__mRun;
            
        for(int j=0; j<PFJetsCHS__; j++)
        {
            
            p4pf.SetPxPyPzE(PFJetsCHS__P4__fCoordinates_fX[j],PFJetsCHS__P4__fCoordinates_fY[j],
                            PFJetsCHS__P4__fCoordinates_fZ[j],PFJetsCHS__P4__fCoordinates_fT[j]);
      
            ////////////Calculation of correction ///
            
            pf_jtpt[j]   = p4pf.Pt();
            pf_jt_eta[j] = p4pf.Eta();
            
            float rho=EvtHdr__mPFRho;
            double area=PFJetsCHS__area_[j];
            double ptraw = pf_jtpt[j]/PFJetsCHS__cor_[j]; //undo islemi
            
            //double ptraw = pf_jtpt[j]+rho*area; //without undo
            pfchs_jec->setJetEta(pf_jt_eta[j]);
            pfchs_jec->setJetPt(ptraw);
            pfchs_jec->setRho(rho);
            pfchs_jec->setJetA(area);
            
            pfchs_ptcorr[j]= ptraw*pfchs_jec->getCorrection();
            
            double pf_pt=pfchs_ptcorr[j];
            double pf_jteta=pf_jt_eta[j];
            
            ///////=====================the end of correction////////
            
      
            
            if( PFJetsCHS__tightID_[j] && (PFMet__et_<0.3*PFMet__sumEt_) && (PFJetsCHS__nhf_[j]<0.9) && (PFJetsCHS__nemf_[j]<0.9) &&(PFJetsCHS__muf_[j]<0.9))
            {
                for ( int k = 0; k<netabins; ++k)
                {
                    
                    if (fabs(pf_jteta)>= etabins[k] && fabs(pf_jteta)< etabins[k+1])
                    {
                        if (k>=0 && k<2){vecttrigcuts=vecttrigcuts1;}//Barrel_trigger_etas
                        if (k>=2 && k<6){vecttrigcuts=vecttrigcuts2;}//EndCap_trigger_etas
                        if (k>=6 && k<8){vecttrigcuts=vecttrigcuts3;}//HF_trigger_etas
                        
                        for(int trnameindex=1; trnameindex<4 ; trnameindex++)
                        {
                            /*if (trnameindex==1 ){ runnocount2->Fill(dummy);
                                if(counter==0){counter=counter+1;
                                //runfile<< EvtHdr__mRun<<" "<<jentry <<endl; runnocount->Fill(dummy);}
                            }*/
                            for(int trdecindex=0; trdecindex<TriggerDecision_.size(); trdecindex++)
                            {
                                
                                if(TriggerDecision_[trdecindex]==trnameindex && ((pf_pt)>=vecttrigcuts[trnameindex]&& (pf_pt)<vecttrigcuts[trnameindex+1]))
                                {
                                    
                                    double w=vectlumi[trnameindex];
                                    //vectptt[trnameindex][k]->Fill(pf_pt);
                                    vectptt[trnameindex][k]->Fill(pf_pt, 1./w);
                                }
                            } //trdecindex
                        } //trnameindex
                    } // fabs
                } //etabin
            } //tightid
        } //jetler
        // counter=0;
    } //event
    myFile.cd();
    myFile.Write();
    myFile.Close();
    
    //time (&end); dif = difftime (end,start); cout<< endl<< "zaman:"<< dif<<endl;
    //runfile.close();
    
}


