//#include "LEGlowpu.C"
void Run()
{
    gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
    gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
    gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
    gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");
    gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
    gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");
    //gROOT->ProcessLine(".exception");
    
    gROOT->ProcessLine(".L newLEGlowpu.C++g");
    gROOT->ProcessLine("newLEGlowpu t");
    //gROOT->ProcessLine(".L ak8lowpu.C++g");
    //gROOT->ProcessLine("ak8lowpu t");
    gROOT->ProcessLine("t.Loop()");
    
    
}
