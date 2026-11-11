#define HSCPSelector_cxx

#include "HSCPSelector.h"
#include <TH2.h>
#include <TH1.h>
#include <functional>
#include <cmath>
#include <cstdlib>
#include <map>
#include <sstream>
#include <stdexcept>


//K and C values are set based on the sample name (see SlaveBegin)
float K(2.54), C(3.14); 
float K_bckg2024(2.83894), C_bckg2024(3.01756); //MC background 2024
float K_data2024(2.8202), C_data2024(2.9784); //Data 2024

// glu pion:
//float K_bckg2024(2.07844), C_bckg2024(2.9784);
//float K_data2024(2.07844), C_data2024(2.9784);


//ADD-SELECTION-METHODS
bool HSCPSelector::PassHSCPpresel_METanalysis_TestPUppiMETCut_Eta2p4(int i){
   if (i<0 || i>=(int)Pt.GetSize()) {
      cout << i << endl;
      return false;
   }
   return (( (Flag_allMETFilters[0] == true) && (RecoPUppiMET[0] > 150.) && (Pt[i] > 50.0) && (Pt_pseudo[i] > 50.0) && (std::abs(Eta[i]) < 2.4) && (NbPixelHit_noL1[i] >= 2) && (FracOfValidHit[i] > 0.8) && (NOM_noL1[i] >= 10) && (isHighPurityTrack[i] == true) && (normChi2[i] < 5.0) && (std::abs(dz[i]) < 0.1) && (std::abs(dxy[i]) < 0.02) && (miniRelIsoAll[i] < 0.02) && (EoP[i] < 0.3) && (IsoSumPt_dr03[i] < 15) && (ptOverptErrptErr[i] < 0.0008) && (Fpix[i] > 0.3) && (ptOverptErr[i] < 1) && (Ih_Strip[i] > C) ));
}

bool HSCPSelector::PassHSCPpresel_METanalysis_TestPUppiMETCut_Eta1(int i){
   if (i<0 || i>=(int)Pt.GetSize()) {
      cout << i << endl;
      return false;
   }
   return (( (Flag_allMETFilters[0] == true) && (RecoPUppiMET[0] > 150.) && (Pt[i] > 50.0) && (Pt_pseudo[i] > 50.0) && (std::abs(Eta[i]) < 1) && (NbPixelHit_noL1[i] >= 2) && (FracOfValidHit[i] > 0.8) && (NOM_noL1[i] >= 10) && (isHighPurityTrack[i] == true) && (normChi2[i] < 5.0) && (std::abs(dz[i]) < 0.1) && (std::abs(dxy[i]) < 0.02) && (miniRelIsoAll[i] < 0.02) && (EoP[i] < 0.3) && (IsoSumPt_dr03[i] < 15) && (ptOverptErrptErr[i] < 0.0008) && (Fpix[i] > 0.3) && (ptOverptErr[i] < 1) && (Ih_Strip[i] > C) ));
}

bool HSCPSelector::PassHSCPpresel_METanalysis_TestPUppiMETCut_Eta1_2p4(int i){
   if (i<0 || i>=(int)Pt.GetSize()) {
      cout << i << endl;
      return false;
   }
   return (( (Flag_allMETFilters[0] == true) && (RecoPUppiMET[0] > 150.) && (Pt[i] > 50.0) && (Pt_pseudo[i] > 50.0) && (std::abs(Eta[i]) >= 1 && std::abs(Eta[i]) < 2.4) && (NbPixelHit_noL1[i] >= 2) && (FracOfValidHit[i] > 0.8) && (NOM_noL1[i] >= 10) && (isHighPurityTrack[i] == true) && (normChi2[i] < 5.0) && (std::abs(dz[i]) < 0.1) && (std::abs(dxy[i]) < 0.02) && (miniRelIsoAll[i] < 0.02) && (EoP[i] < 0.3) && (IsoSumPt_dr03[i] < 15) && (ptOverptErrptErr[i] < 0.0008) && (Fpix[i] > 0.3) && (ptOverptErr[i] < 1) && (Ih_Strip[i] > C) ));
}

bool HSCPSelector::PassHSCPpresel_METanalysis_TestPUppiMETCut_EoP_0p1_Eta2p4(int i){
   if (i<0 || i>=(int)Pt.GetSize()) {
      cout << i << endl;
      return false;
   }
   return (( (Flag_allMETFilters[0] == true) && (RecoPUppiMET[0] > 150.) && (Pt[i] > 50.0) && (Pt_pseudo[i] > 50.0) && (std::abs(Eta[i]) < 2.4) && (NbPixelHit_noL1[i] >= 2) && (FracOfValidHit[i] > 0.8) && (NOM_noL1[i] >= 10) && (isHighPurityTrack[i] == true) && (normChi2[i] < 5.0) && (std::abs(dz[i]) < 0.1) && (std::abs(dxy[i]) < 0.02) && (miniRelIsoAll[i] < 0.02) && (EoP[i] < 0.1) && (IsoSumPt_dr03[i] < 15) && (ptOverptErrptErr[i] < 0.0008) && (Fpix[i] > 0.3) && (ptOverptErr[i] < 1) && (Ih_Strip[i] > C) ));
}

bool HSCPSelector::PassHSCPpresel_METanalysis_TestPUppiMETCut_EoP_0p1_Eta1(int i){
   if (i<0 || i>=(int)Pt.GetSize()) {
      cout << i << endl;
      return false;
   }
   return (( (Flag_allMETFilters[0] == true) && (RecoPUppiMET[0] > 150.) && (Pt[i] > 50.0) && (Pt_pseudo[i] > 50.0) && (std::abs(Eta[i]) < 1) && (NbPixelHit_noL1[i] >= 2) && (FracOfValidHit[i] > 0.8) && (NOM_noL1[i] >= 10) && (isHighPurityTrack[i] == true) && (normChi2[i] < 5.0) && (std::abs(dz[i]) < 0.1) && (std::abs(dxy[i]) < 0.02) && (miniRelIsoAll[i] < 0.02) && (EoP[i] < 0.1) && (IsoSumPt_dr03[i] < 15) && (ptOverptErrptErr[i] < 0.0008) && (Fpix[i] > 0.3) && (ptOverptErr[i] < 1) && (Ih_Strip[i] > C) ));
}

bool HSCPSelector::PassHSCPpresel_METanalysis_TestPUppiMETCut_EoP_0p1_Eta1_2p4(int i){
   if (i<0 || i>=(int)Pt.GetSize()) {
      cout << i << endl;
      return false;
   }
   return (( (Flag_allMETFilters[0] == true) && (RecoPUppiMET[0] > 150.) && (Pt[i] > 50.0) && (Pt_pseudo[i] > 50.0) && (std::abs(Eta[i]) >= 1 && std::abs(Eta[i]) < 2.4) && (NbPixelHit_noL1[i] >= 2) && (FracOfValidHit[i] > 0.8) && (NOM_noL1[i] >= 10) && (isHighPurityTrack[i] == true) && (normChi2[i] < 5.0) && (std::abs(dz[i]) < 0.1) && (std::abs(dxy[i]) < 0.02) && (miniRelIsoAll[i] < 0.02) && (EoP[i] < 0.1) && (IsoSumPt_dr03[i] < 15) && (ptOverptErrptErr[i] < 0.0008) && (Fpix[i] > 0.3) && (ptOverptErr[i] < 1) && (Ih_Strip[i] > C) ));
}

bool HSCPSelector::PassHSCPpresel_METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5_Eta2p4(int i){
   if (i<0 || i>=(int)Pt.GetSize()) {
      cout << i << endl;
      return false;
   }
   return (( (Flag_allMETFilters[0] == true) && (RecoPUppiMET[0] > 150.) && (Pt[i] > 50.0) && (Pt_pseudo[i] > 50.0) && (std::abs(Eta[i]) < 2.4) && (NbPixelHit_noL1[i] >= 2) && (FracOfValidHit[i] > 0.8) && (NOM_noL1[i] >= 10) && (isHighPurityTrack[i] == true) && (normChi2[i] < 5.0) && (std::abs(dz[i]) < 0.1) && (std::abs(dxy[i]) < 0.02) && (miniRelIsoAll[i] < 0.02) && (EoP[i] < 0.3) && (IsoSumPt_dr03[i] < 15) && (ptOverptErrptErr[i] < 0.0008) && (Fpix[i] > 0.3) && (ptOverptErr[i] < 0.5) && (Ih_Strip[i] > C) ));
}

bool HSCPSelector::PassHSCPpresel_METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5_Eta1(int i){
   if (i<0 || i>=(int)Pt.GetSize()) {
      cout << i << endl;
      return false;
   }
   return (( (Flag_allMETFilters[0] == true) && (RecoPUppiMET[0] > 150.) && (Pt[i] > 50.0) && (Pt_pseudo[i] > 50.0) && (std::abs(Eta[i]) < 1) && (NbPixelHit_noL1[i] >= 2) && (FracOfValidHit[i] > 0.8) && (NOM_noL1[i] >= 10) && (isHighPurityTrack[i] == true) && (normChi2[i] < 5.0) && (std::abs(dz[i]) < 0.1) && (std::abs(dxy[i]) < 0.02) && (miniRelIsoAll[i] < 0.02) && (EoP[i] < 0.3) && (IsoSumPt_dr03[i] < 15) && (ptOverptErrptErr[i] < 0.0008) && (Fpix[i] > 0.3) && (ptOverptErr[i] < 0.5) && (Ih_Strip[i] > C) ));
}

bool HSCPSelector::PassHSCPpresel_METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5_Eta1_2p4(int i){
   if (i<0 || i>=(int)Pt.GetSize()) {
      cout << i << endl;
      return false;
   }
   return (( (Flag_allMETFilters[0] == true) && (RecoPUppiMET[0] > 150.) && (Pt[i] > 50.0) && (Pt_pseudo[i] > 50.0) && (std::abs(Eta[i]) >= 1 && std::abs(Eta[i]) < 2.4) && (NbPixelHit_noL1[i] >= 2) && (FracOfValidHit[i] > 0.8) && (NOM_noL1[i] >= 10) && (isHighPurityTrack[i] == true) && (normChi2[i] < 5.0) && (std::abs(dz[i]) < 0.1) && (std::abs(dxy[i]) < 0.02) && (miniRelIsoAll[i] < 0.02) && (EoP[i] < 0.3) && (IsoSumPt_dr03[i] < 15) && (ptOverptErrptErr[i] < 0.0008) && (Fpix[i] > 0.3) && (ptOverptErr[i] < 0.5) && (Ih_Strip[i] > C) ));
}

bool HSCPSelector::PassHSCPpresel_METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5_EoP_0p1_Eta2p4(int i){
   if (i<0 || i>=(int)Pt.GetSize()) {
      cout << i << endl;
      return false;
   }
   return (( (Flag_allMETFilters[0] == true) && (RecoPUppiMET[0] > 150.) && (Pt[i] > 50.0) && (Pt_pseudo[i] > 50.0) && (std::abs(Eta[i]) < 2.4) && (NbPixelHit_noL1[i] >= 2) && (FracOfValidHit[i] > 0.8) && (NOM_noL1[i] >= 10) && (isHighPurityTrack[i] == true) && (normChi2[i] < 5.0) && (std::abs(dz[i]) < 0.1) && (std::abs(dxy[i]) < 0.02) && (miniRelIsoAll[i] < 0.02) && (EoP[i] < 0.1) && (IsoSumPt_dr03[i] < 15) && (ptOverptErrptErr[i] < 0.0008) && (Fpix[i] > 0.3) && (ptOverptErr[i] < 0.5) && (Ih_Strip[i] > C) ));
}

bool HSCPSelector::PassHSCPpresel_METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5_EoP_0p1_Eta1(int i){
   if (i<0 || i>=(int)Pt.GetSize()) {
      cout << i << endl;
      return false;
   }
   return (( (Flag_allMETFilters[0] == true) && (RecoPUppiMET[0] > 150.) && (Pt[i] > 50.0) && (Pt_pseudo[i] > 50.0) && (std::abs(Eta[i]) < 1) && (NbPixelHit_noL1[i] >= 2) && (FracOfValidHit[i] > 0.8) && (NOM_noL1[i] >= 10) && (isHighPurityTrack[i] == true) && (normChi2[i] < 5.0) && (std::abs(dz[i]) < 0.1) && (std::abs(dxy[i]) < 0.02) && (miniRelIsoAll[i] < 0.02) && (EoP[i] < 0.1) && (IsoSumPt_dr03[i] < 15) && (ptOverptErrptErr[i] < 0.0008) && (Fpix[i] > 0.3) && (ptOverptErr[i] < 0.5) && (Ih_Strip[i] > C) ));
}

bool HSCPSelector::PassHSCPpresel_METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5_EoP_0p1_Eta1_2p4(int i){
   if (i<0 || i>=(int)Pt.GetSize()) {
      cout << i << endl;
      return false;
   }
   return (( (Flag_allMETFilters[0] == true) && (RecoPUppiMET[0] > 150.) && (Pt[i] > 50.0) && (Pt_pseudo[i] > 50.0) && (std::abs(Eta[i]) >= 1 && std::abs(Eta[i]) < 2.4) && (NbPixelHit_noL1[i] >= 2) && (FracOfValidHit[i] > 0.8) && (NOM_noL1[i] >= 10) && (isHighPurityTrack[i] == true) && (normChi2[i] < 5.0) && (std::abs(dz[i]) < 0.1) && (std::abs(dxy[i]) < 0.02) && (miniRelIsoAll[i] < 0.02) && (EoP[i] < 0.1) && (IsoSumPt_dr03[i] < 15) && (ptOverptErrptErr[i] < 0.0008) && (Fpix[i] > 0.3) && (ptOverptErr[i] < 0.5) && (Ih_Strip[i] > C) ));
}



// The option string given by macro.cc is a list of key=value pairs separated by commas:
//   ptcut, etabins, ihbins, pbins, massbins, fpixbins, sample, out (output file), nosel (0 or 1)
std::map<std::string, std::string> HSCPSelector::parseOptions() const {
    std::map<std::string, std::string> opt;
    std::stringstream ss(GetOption());
    for (std::string item; std::getline(ss, item, ','); ) {
        const size_t eq = item.find('=');
        if (eq == std::string::npos) continue;
        opt[item.substr(0, eq)] = item.substr(eq + 1);
    }
    for (const char* key : {"ptcut", "etabins", "ihbins", "pbins", "massbins", "fpixbins", "sample"}) {
        if (!opt.count(key)) throw std::runtime_error(std::string("HSCPSelector: option '") + key + "' missing in \"" + GetOption() + "\"");
    }
    return opt;
}


void HSCPSelector::Begin(TTree *tree) {
    std::map<std::string, std::string> opt = parseOptions();

    oFile_ = (opt.count("out") && !opt["out"].empty()) ? opt["out"] : opt["sample"] + ".root";
    fillNosel_ = !(opt.count("nosel") && opt["nosel"] == "0");


    // Trigger scale-factor tables: directory HSCP_SF_DIR, default ../inc (relative to macros/)
    std::string sfDir = "../inc";
    if (const char* env = std::getenv("HSCP_SF_DIR")) { if (*env) sfDir = env; }
    loadSF2D(sfDir + "/SF_orMETtrg_PUppiMET_VS_PseudoMET__TriggerEffCalib_table_plain.txt", SF2D_PseudoMETlovalue, SF2D_PseudoMEThivalue, SF2D_PUppiMETlovalue, SF2D_PUppiMEThivalue, SF2D_Down, SF2D, SF2D_Up);
    loadSF1D(sfDir + "/SF_PseudoMET.txt", SF1D_PseudoMETvalue, SF1Dpseudo_Down, SF1Dpseudo, SF1Dpseudo_Up);
    loadSF1D(sfDir + "/SF_PUppiMET.txt", SF1D_PUppiMETvalue, SF1Dpuppi_Down, SF1Dpuppi, SF1Dpuppi_Up);


    //FILL-SELECTION-VECTOR
selections_.push_back(&HSCPSelector::PassHSCPpresel_METanalysis_TestPUppiMETCut_Eta2p4);
selLabels_.push_back("METanalysis_TestPUppiMETCut_Eta2p4");

selections_.push_back(&HSCPSelector::PassHSCPpresel_METanalysis_TestPUppiMETCut_Eta1);
selLabels_.push_back("METanalysis_TestPUppiMETCut_Eta1");

selections_.push_back(&HSCPSelector::PassHSCPpresel_METanalysis_TestPUppiMETCut_Eta1_2p4);
selLabels_.push_back("METanalysis_TestPUppiMETCut_Eta1_2p4");

selections_.push_back(&HSCPSelector::PassHSCPpresel_METanalysis_TestPUppiMETCut_EoP_0p1_Eta2p4);
selLabels_.push_back("METanalysis_TestPUppiMETCut_EoP_0p1_Eta2p4");

selections_.push_back(&HSCPSelector::PassHSCPpresel_METanalysis_TestPUppiMETCut_EoP_0p1_Eta1);
selLabels_.push_back("METanalysis_TestPUppiMETCut_EoP_0p1_Eta1");

selections_.push_back(&HSCPSelector::PassHSCPpresel_METanalysis_TestPUppiMETCut_EoP_0p1_Eta1_2p4);
selLabels_.push_back("METanalysis_TestPUppiMETCut_EoP_0p1_Eta1_2p4");

selections_.push_back(&HSCPSelector::PassHSCPpresel_METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5_Eta2p4);
selLabels_.push_back("METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5_Eta2p4");

selections_.push_back(&HSCPSelector::PassHSCPpresel_METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5_Eta1);
selLabels_.push_back("METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5_Eta1");

selections_.push_back(&HSCPSelector::PassHSCPpresel_METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5_Eta1_2p4);
selLabels_.push_back("METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5_Eta1_2p4");

selections_.push_back(&HSCPSelector::PassHSCPpresel_METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5_EoP_0p1_Eta2p4);
selLabels_.push_back("METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5_EoP_0p1_Eta2p4");

selections_.push_back(&HSCPSelector::PassHSCPpresel_METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5_EoP_0p1_Eta1);
selLabels_.push_back("METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5_EoP_0p1_Eta1");

selections_.push_back(&HSCPSelector::PassHSCPpresel_METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5_EoP_0p1_Eta1_2p4);
selLabels_.push_back("METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5_EoP_0p1_Eta1_2p4");


    std::cout << std::endl;
    std::cout << "Number of selections: " << selLabels_.size() << std::endl;
    std::cout << "Output file: " << oFile_.c_str() << std::endl;
    std::cout << std::endl;
}



void HSCPSelector::SlaveBegin(TTree *tree) {
    std::map<std::string, std::string> opt = parseOptions();

    ptcut_    = std::stod(opt["ptcut"]);
    etabins_  = std::stoi(opt["etabins"]);
    ihbins_   = std::stoi(opt["ihbins"]);
    pbins_    = std::stoi(opt["pbins"]);
    massbins_ = std::stoi(opt["massbins"]);
    fpixbins_ = std::stoi(opt["fpixbins"]);
    sample_   = opt["sample"];


    if(sample_.find("2024") != std::string::npos) {
        K = K_data2024;
        C = C_data2024;
    }
    if (sample_.find("Gluino") != std::string::npos ||
        sample_.find("Stau") != std::string::npos ||
        sample_.find("Stop") != std::string::npos ||
        sample_.find("QCD") != std::string::npos || 
        sample_.find("Wjets") != std::string::npos ||
        sample_.find("WjetMuNu") != std::string::npos ||
        sample_.find("TTbar") != std::string::npos) {
        K = K_bckg2024;
        C = C_bckg2024;
    }

    cout << "Value of K: " << K << endl;
    cout << "         C: " << C << endl;


    //-------------------------------------
    //create RegionMassPlot for all selection
    //-------------------------------------
    for(unsigned int i=0; i<selLabels_.size();i++)
    {
        cout << "Initialisation of the regions for selection: " << selLabels_[i] << endl;

        //Create names for each plot in each regions (slices in Fpixels)
        
        std::string label_FpixA_3f4 = regFpixA_3f4 + "_" + selLabels_[i]; 
        std::string label_FpixA_3f6 = regFpixA_3f6 + "_" + selLabels_[i]; 
        std::string label_FpixA_3f8 = regFpixA_3f8 + "_" + selLabels_[i]; 
        std::string label_FpixA_3f9 = regFpixA_3f9 + "_" + selLabels_[i]; 
        std::string label_FpixA_4f5 = regFpixA_4f5 + "_" + selLabels_[i]; 
        std::string label_FpixA_5f6 = regFpixA_5f6 + "_" + selLabels_[i]; 
        std::string label_FpixA_6f7 = regFpixA_6f7 + "_" + selLabels_[i]; 
        std::string label_FpixA_6f9 = regFpixA_6f9 + "_" + selLabels_[i];
        std::string label_FpixA_7f8 = regFpixA_7f8 + "_" + selLabels_[i]; 
        std::string label_FpixA_8f9 = regFpixA_8f9 + "_" + selLabels_[i]; 
        std::string label_FpixA_9f10 = regFpixA_9f10 + "_" + selLabels_[i];
        std::string label_FpixA_99f10 = regFpixA_99f10 + "_" + selLabels_[i];
        std::string label_FpixA_999f10 = regFpixA_999f10 + "_" + selLabels_[i];

        std::string label_FpixB_3f4 = regFpixB_3f4 + "_" + selLabels_[i]; 
        std::string label_FpixB_3f6 = regFpixB_3f6 + "_" + selLabels_[i]; 
        std::string label_FpixB_3f8 = regFpixB_3f8 + "_" + selLabels_[i]; 
        std::string label_FpixB_3f9 = regFpixB_3f9 + "_" + selLabels_[i]; 
        std::string label_FpixB_4f5 = regFpixB_4f5 + "_" + selLabels_[i]; 
        std::string label_FpixB_5f6 = regFpixB_5f6 + "_" + selLabels_[i]; 
        std::string label_FpixB_6f7 = regFpixB_6f7 + "_" + selLabels_[i]; 
        std::string label_FpixB_6f9 = regFpixB_6f9 + "_" + selLabels_[i]; 
        std::string label_FpixB_7f8 = regFpixB_7f8 + "_" + selLabels_[i]; 
        std::string label_FpixB_8f9 = regFpixB_8f9 + "_" + selLabels_[i]; 
        std::string label_FpixB_8f10 = regFpixB_8f10 + "_" + selLabels_[i];
        std::string label_FpixB_9f10 = regFpixB_9f10 + "_" + selLabels_[i];
        std::string label_FpixB_99f10 = regFpixB_99f10 + "_" + selLabels_[i];
        std::string label_FpixB_999f10 = regFpixB_999f10 + "_" + selLabels_[i];

        std::string label_FpixC_3f4 = regFpixC_3f4 + "_" + selLabels_[i]; 
        std::string label_FpixC_3f6 = regFpixC_3f6 + "_" + selLabels_[i]; 
        std::string label_FpixC_3f8 = regFpixC_3f8 + "_" + selLabels_[i]; 
        std::string label_FpixC_3f9 = regFpixC_3f9 + "_" + selLabels_[i]; 
        std::string label_FpixC_4f5 = regFpixC_4f5 + "_" + selLabels_[i]; 
        std::string label_FpixC_5f6 = regFpixC_5f6 + "_" + selLabels_[i]; 
        std::string label_FpixC_6f7 = regFpixC_6f7 + "_" + selLabels_[i]; 
        std::string label_FpixC_6f9 = regFpixC_6f9 + "_" + selLabels_[i]; 
        std::string label_FpixC_7f8 = regFpixC_7f8 + "_" + selLabels_[i]; 
        std::string label_FpixC_8f9 = regFpixC_8f9 + "_" + selLabels_[i]; 

        std::string label_FpixD_3f4 = regFpixD_3f4 + "_" + selLabels_[i];
        std::string label_FpixD_3f8 = regFpixD_3f8 + "_" + selLabels_[i];
        std::string label_FpixD_4f5 = regFpixD_4f5 + "_" + selLabels_[i]; 
        std::string label_FpixD_5f6 = regFpixD_5f6 + "_" + selLabels_[i]; 
        std::string label_FpixD_6f7 = regFpixD_6f7 + "_" + selLabels_[i]; 
        std::string label_FpixD_6f9 = regFpixD_6f9 + "_" + selLabels_[i];
        std::string label_FpixD_7f8 = regFpixD_7f8 + "_" + selLabels_[i]; 
        std::string label_FpixD_8f9 = regFpixD_8f9 + "_" + selLabels_[i]; 
        std::string label_FpixD_8f10 = regFpixD_8f10 + "_" + selLabels_[i];
        std::string label_FpixD_9f10 = regFpixD_9f10 + "_" + selLabels_[i]; 
        std::string label_FpixD_99f10 = regFpixD_99f10 + "_" + selLabels_[i]; 
        std::string label_FpixD_999f10 = regFpixD_999f10 + "_" + selLabels_[i]; 


        //Create objects RegionMassPlot using the names defined above

        RegionMassPlot regA_3f4(label_FpixA_3f4.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C);
        RegionMassPlot regA_3f6(label_FpixA_3f6.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regA_3f8(label_FpixA_3f8.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regA_3f9(label_FpixA_3f9.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regA_4f5(label_FpixA_4f5.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regA_5f6(label_FpixA_5f6.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regA_6f7(label_FpixA_6f7.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C);
        RegionMassPlot regA_6f9(label_FpixA_6f9.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C);
        RegionMassPlot regA_7f8(label_FpixA_7f8.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C);
        RegionMassPlot regA_8f9(label_FpixA_8f9.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regA_9f10(label_FpixA_9f10.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regA_99f10(label_FpixA_99f10.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regA_999f10(label_FpixA_999f10.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 

        RegionMassPlot regB_3f4(label_FpixB_3f4.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regB_3f6(label_FpixB_3f6.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regB_3f8(label_FpixB_3f8.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regB_3f9(label_FpixB_3f9.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regB_4f5(label_FpixB_4f5.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regB_5f6(label_FpixB_5f6.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regB_6f7(label_FpixB_6f7.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regB_6f9(label_FpixB_6f9.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regB_7f8(label_FpixB_7f8.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C);     
        RegionMassPlot regB_8f9(label_FpixB_8f9.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regB_8f10(label_FpixB_8f10.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regB_9f10(label_FpixB_9f10.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regB_99f10(label_FpixB_99f10.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regB_999f10(label_FpixB_999f10.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C);

        RegionMassPlot regC_3f4(label_FpixC_3f4.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regC_3f6(label_FpixC_3f6.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regC_3f8(label_FpixC_3f8.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regC_3f9(label_FpixC_3f9.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regC_4f5(label_FpixC_4f5.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regC_5f6(label_FpixC_5f6.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regC_6f7(label_FpixC_6f7.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regC_6f9(label_FpixC_6f9.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regC_7f8(label_FpixC_7f8.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regC_8f9(label_FpixC_8f9.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 

        RegionMassPlot regD_3f4(label_FpixD_3f4.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regD_3f8(label_FpixD_3f8.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regD_4f5(label_FpixD_4f5.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regD_5f6(label_FpixD_5f6.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regD_6f7(label_FpixD_6f7.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regD_6f9(label_FpixD_6f9.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regD_7f8(label_FpixD_7f8.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C);     
        RegionMassPlot regD_8f9(label_FpixD_8f9.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regD_8f10(label_FpixD_8f10.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regD_9f10(label_FpixD_9f10.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regD_99f10(label_FpixD_99f10.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
        RegionMassPlot regD_999f10(label_FpixD_999f10.c_str(), etabins_, ihbins_, pbins_, massbins_, fpixbins_, C); 
    

        //push the RegionMassPlot objects created in vectors, because you have a set of different plot for each region, for each selection (defined in the configuration file) 
        
        vmrp_regionA_3f4.push_back(std::move(regA_3f4));
        vmrp_regionA_3f6.push_back(std::move(regA_3f6));
        vmrp_regionA_3f8.push_back(std::move(regA_3f8));
        vmrp_regionA_3f9.push_back(std::move(regA_3f9));
        vmrp_regionA_4f5.push_back(std::move(regA_4f5));
        vmrp_regionA_5f6.push_back(std::move(regA_5f6));
        vmrp_regionA_6f7.push_back(std::move(regA_6f7));
        vmrp_regionA_6f9.push_back(std::move(regA_6f9));
        vmrp_regionA_7f8.push_back(std::move(regA_7f8));
        vmrp_regionA_8f9.push_back(std::move(regA_8f9));
        vmrp_regionA_9f10.push_back(std::move(regA_9f10));
        vmrp_regionA_99f10.push_back(std::move(regA_99f10));
        vmrp_regionA_999f10.push_back(std::move(regA_999f10));

        vmrp_regionB_3f4.push_back(std::move(regB_3f4));
        vmrp_regionB_3f6.push_back(std::move(regB_3f6));
        vmrp_regionB_3f8.push_back(std::move(regB_3f8));
        vmrp_regionB_3f9.push_back(std::move(regB_3f9));
        vmrp_regionB_4f5.push_back(std::move(regB_4f5));
        vmrp_regionB_5f6.push_back(std::move(regB_5f6));
        vmrp_regionB_6f7.push_back(std::move(regB_6f7));
        vmrp_regionB_6f9.push_back(std::move(regB_6f9));
        vmrp_regionB_7f8.push_back(std::move(regB_7f8));
        vmrp_regionB_8f9.push_back(std::move(regB_8f9));
        vmrp_regionB_8f10.push_back(std::move(regB_8f10));
        vmrp_regionB_9f10.push_back(std::move(regB_9f10));
        vmrp_regionB_99f10.push_back(std::move(regB_99f10));
        vmrp_regionB_999f10.push_back(std::move(regB_999f10));

        vmrp_regionC_3f4.push_back(std::move(regC_3f4));
        vmrp_regionC_3f6.push_back(std::move(regC_3f6));
        vmrp_regionC_3f8.push_back(std::move(regC_3f8));
        vmrp_regionC_3f9.push_back(std::move(regC_3f9));    
        vmrp_regionC_4f5.push_back(std::move(regC_4f5));
        vmrp_regionC_5f6.push_back(std::move(regC_5f6));
        vmrp_regionC_6f7.push_back(std::move(regC_6f7));
        vmrp_regionC_6f9.push_back(std::move(regC_6f9));
        vmrp_regionC_7f8.push_back(std::move(regC_7f8));
        vmrp_regionC_8f9.push_back(std::move(regC_8f9));

        vmrp_regionD_3f4.push_back(std::move(regD_3f4));
        vmrp_regionD_3f8.push_back(std::move(regD_3f8));
        vmrp_regionD_4f5.push_back(std::move(regD_4f5));
        vmrp_regionD_5f6.push_back(std::move(regD_5f6));
        vmrp_regionD_6f7.push_back(std::move(regD_6f7));
        vmrp_regionD_6f9.push_back(std::move(regD_6f9));
        vmrp_regionD_7f8.push_back(std::move(regD_7f8));
        vmrp_regionD_8f9.push_back(std::move(regD_8f9));
        vmrp_regionD_8f10.push_back(std::move(regD_8f10));
        vmrp_regionD_9f10.push_back(std::move(regD_9f10));
        vmrp_regionD_99f10.push_back(std::move(regD_99f10));
        vmrp_regionD_999f10.push_back(std::move(regD_999f10));


        CPlots plots;

        // mass calibration
        // plots.AddHisto2F(selLabels_[i] + "_Mu_ih_vs_betagamma", 500, 0., 15000., 100, 0., 10.);
        // plots.AddHisto2F(selLabels_[i] + "_Mu_ih_vs_p", 500, 0., 6000., 100, 0.,10.);
        // plots.AddHisto2F(selLabels_[i] + "_KPpi_ih_vs_betagamma", 500, 0., 5000., 100, 0., 10.);
        // plots.AddHisto2F(selLabels_[i] + "_KPpi_ih_vs_p", 500, 0., 3000., 100, 0., 10.);
        // plots.AddHisto2F(selLabels_[i] + "_Gluino_ih_vs_betagamma", 1000, 0., 5., 500, 0., 50.);
        // plots.AddHisto2F(selLabels_[i] + "_Gluino_ih_vs_p", 500, 0., 3000., 500, 0.,50.);


        // Trigger efficiency
        plots.AddHisto1F(selLabels_[i] + "_L1MET", 40, 0, 1000);
        plots.AddHisto1F(selLabels_[i] + "_HLTCaloMET", 40, 0, 1000);
        plots.AddHisto1F(selLabels_[i] + "_HLTCaloMHT", 40, 0, 1000);
        plots.AddHisto1F(selLabels_[i] + "_HLTPFMHT", 40, 0, 1000);
        plots.AddHisto1F(selLabels_[i] + "_HLTPFMET", 40, 0, 1000);


        plots.AddHisto1F(selLabels_[i] + "_PFtrackPT", 300, 0, 3000);
        plots.AddHisto1F(selLabels_[i] + "_if___HLT_MET105_IsoTrk50___PFtrackPT", 300, 0, 3000);
        
        plots.AddHisto2F(selLabels_[i] + "_PUppiMET_VS_PseudoMET", 100, 0, 2500, 100, 0, 2500);
        plots.AddHisto2F(selLabels_[i] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMET_VS_PseudoMET", 100, 0, 2500, 100, 0, 2500);
        plots.AddHisto2F(selLabels_[i] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMET_VS_PseudoMET", 100, 0, 2500, 100, 0, 2500);
        plots.AddHisto2F(selLabels_[i] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMET_VS_PseudoMET", 100, 0, 2500, 100, 0, 2500);
        plots.AddHisto2F(selLabels_[i] + "_if___HLT_MET105_IsoTrk50___PUppiMET_VS_PseudoMET", 100, 0, 2500, 100, 0, 2500);
        plots.AddHisto2F(selLabels_[i] + "_if___orMETtrg___PUppiMET_VS_PseudoMET", 100, 0, 2500, 100, 0, 2500);
        plots.AddHisto2F(selLabels_[i] + "_if___orMET3a4trg1___PUppiMET_VS_PseudoMET", 100, 0, 2500, 100, 0, 2500);
        plots.AddHisto2F(selLabels_[i] + "_if___orMET3a4trg2___PUppiMET_VS_PseudoMET", 100, 0, 2500, 100, 0, 2500);
        plots.AddHisto2F(selLabels_[i] + "_if___orMET3a4trg3___PUppiMET_VS_PseudoMET", 100, 0, 2500, 100, 0, 2500);
        plots.AddHisto2F(selLabels_[i] + "_if___orMET3a4trg4___PUppiMET_VS_PseudoMET", 100, 0, 2500, 100, 0, 2500);

        plots.AddHisto2F(selLabels_[i] + "_PUppiMETNoMu_VS_PseudoMET", 100, 0, 2500, 100, 0, 2500);
        plots.AddHisto2F(selLabels_[i] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMETNoMu_VS_PseudoMET", 100, 0, 2500, 100, 0, 2500);
        plots.AddHisto2F(selLabels_[i] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMETNoMu_VS_PseudoMET", 100, 0, 2500, 100, 0, 2500);
        plots.AddHisto2F(selLabels_[i] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMETNoMu_VS_PseudoMET", 100, 0, 2500, 100, 0, 2500);
        plots.AddHisto2F(selLabels_[i] + "_if___HLT_MET105_IsoTrk50___PUppiMETNoMu_VS_PseudoMET", 100, 0, 2500, 100, 0, 2500);
        plots.AddHisto2F(selLabels_[i] + "_if___orMETtrg___PUppiMETNoMu_VS_PseudoMET", 100, 0, 2500, 100, 0, 2500);
        plots.AddHisto2F(selLabels_[i] + "_if___orMET3a4trg1___PUppiMETNoMu_VS_PseudoMET", 100, 0, 2500, 100, 0, 2500);
        plots.AddHisto2F(selLabels_[i] + "_if___orMET3a4trg2___PUppiMETNoMu_VS_PseudoMET", 100, 0, 2500, 100, 0, 2500);
        plots.AddHisto2F(selLabels_[i] + "_if___orMET3a4trg3___PUppiMETNoMu_VS_PseudoMET", 100, 0, 2500, 100, 0, 2500);
        plots.AddHisto2F(selLabels_[i] + "_if___orMET3a4trg4___PUppiMETNoMu_VS_PseudoMET", 100, 0, 2500, 100, 0, 2500);


        plots.AddHisto1F(selLabels_[i] + "_PUppiMET__SF", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_PUppiMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_PUppiMET__postTrigger__SF", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_PUppiMET__postTrigger", 100, 0, 2500);

        plots.AddHisto1F(selLabels_[i] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___HLT_MET105_IsoTrk50___PUppiMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___orMETtrg___PUppiMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___orMET3a4trg1___PUppiMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___orMET3a4trg2___PUppiMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___orMET3a4trg3___PUppiMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___orMET3a4trg4___PUppiMET", 100, 0, 2500);

        plots.AddHisto1F(selLabels_[i] + "_PUppiMETNoMu", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMETNoMu", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMETNoMu", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMETNoMu", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___HLT_MET105_IsoTrk50___PUppiMETNoMu", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___orMETtrg___PUppiMETNoMu", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___orMET3a4trg1___PUppiMETNoMu", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___orMET3a4trg2___PUppiMETNoMu", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___orMET3a4trg3___PUppiMETNoMu", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___orMET3a4trg4___PUppiMETNoMu", 100, 0, 2500);
        
        plots.AddHisto1F(selLabels_[i] + "_PseudoCaloMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___HLT_MET105_IsoTrk50___PseudoCaloMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___orMETtrg___PseudoCaloMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___orMET3a4trg1___PseudoCaloMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___orMET3a4trg2___PseudoCaloMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___orMET3a4trg3___PseudoCaloMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___orMET3a4trg4___PseudoCaloMET", 100, 0, 2500);

        plots.AddHisto1F(selLabels_[i] + "_RecoPFMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___HLT_MET105_IsoTrk50___RecoPFMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___orMETtrg___RecoPFMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___orMET3a4trg1___RecoPFMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___orMET3a4trg2___RecoPFMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___orMET3a4trg3___RecoPFMET", 100, 0, 2500);
        plots.AddHisto1F(selLabels_[i] + "_if___orMET3a4trg4___RecoPFMET", 100, 0, 2500);

        plots.AddHisto1F(selLabels_[i] + "_PseudoCaloMET__RecoPFMETCut", 100, 0, 2500); 
        plots.AddHisto1F(selLabels_[i] + "_if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET__RecoPFMETCut", 100, 0, 2500); 
        plots.AddHisto1F(selLabels_[i] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET__RecoPFMETCut", 100, 0, 2500); 
        plots.AddHisto1F(selLabels_[i] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET__RecoPFMETCut", 100, 0, 2500); 
        plots.AddHisto1F(selLabels_[i] + "_if___HLT_MET105_IsoTrk50___PseudoCaloMET__RecoPFMETCut", 100, 0, 2500); 
        plots.AddHisto1F(selLabels_[i] + "_if___orMETtrg___PseudoCaloMET__RecoPFMETCut", 100, 0, 2500); 

        plots.AddHisto1F(selLabels_[i] + "_RecoPFMET__PseudoCaloMETCut", 100, 0, 2500); 
        plots.AddHisto1F(selLabels_[i] + "_if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET__PseudoCaloMETCut", 100, 0, 2500); 
        plots.AddHisto1F(selLabels_[i] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET__PseudoCaloMETCut", 100, 0, 2500);         
        plots.AddHisto1F(selLabels_[i] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET__PseudoCaloMETCut", 100, 0, 2500); 
        plots.AddHisto1F(selLabels_[i] + "_if___HLT_MET105_IsoTrk50___RecoPFMET__PseudoCaloMETCut", 100, 0, 2500); 
        plots.AddHisto1F(selLabels_[i] + "_if___orMETtrg___RecoPFMET__PseudoCaloMETCut", 100, 0, 2500); 
        

        // General plots
        plots.AddHisto1F(selLabels_[i] + "_ndEdx_StripOnly", 50, 0, 50);
        plots.AddHisto1F(selLabels_[i] + "_Ptpseudo", 300, 0, 3000);
        plots.AddHisto1F(selLabels_[i] + "_P", 200, 0, 8000);
        plots.AddHisto1F(selLabels_[i] + "_P2", 200, 0, 8000);
        plots.AddHisto1F(selLabels_[i] + "_Psfup", 200, 0, 8000);
        plots.AddHisto1F(selLabels_[i] + "_Psfdown", 200, 0, 8000);
        plots.AddHisto1F(selLabels_[i] + "_pT", 300, 0, 3000);
        plots.AddHisto1F(selLabels_[i] + "_10000oP", 300, 0, 300);
        plots.AddHisto2F(selLabels_[i] + "_10000oP_vs_Eta", 300, 0, 300, 60, -3, +3);
        plots.AddHisto1F(selLabels_[i] + "_eta", 60, -3, +3);
        plots.AddHisto1F(selLabels_[i] + "_phi", 64, -3.2, 3.2);
        plots.AddHisto1F(selLabels_[i] + "_PFMiniIso", 1000, 0, 0.5);
        plots.AddHisto1F(selLabels_[i] + "_TrkIso", 400, 0, 200);
        plots.AddHisto1F(selLabels_[i] + "_EoverP", 500, 0, 5);
        plots.AddHisto1F(selLabels_[i] + "_PtErr_over_PtPt", 500, 0, 0.01);
        plots.AddHisto1F(selLabels_[i] + "_Fpix", 21, 0, 1.1);
        plots.AddHisto1F(selLabels_[i] + "_PtErr_over_Pt", 800, 0, 4);
        plots.AddHisto1F(selLabels_[i] + "_Ih", 600, 0, 30);
        plots.AddHisto1F(selLabels_[i] + "_Ih_rescaled", 200, 0, 10);
        plots.AddHisto1F(selLabels_[i] + "_Ih_oldCorr", 200, 0, 10);
        plots.AddHisto1F(selLabels_[i] + "_Ih_noSF", 200, 0, 10);
        plots.AddHisto1F(selLabels_[i] + "_PthatQCD", 2000, 0, 2000);

        plots.AddHisto2F(selLabels_[i] + "_Fpix_vs_RunNumber", 5000, 378000, 388000, 20, 0, 1);
        plots.AddHisto2F(selLabels_[i] + "_Ih_vs_RunNumber", 5000, 378000, 388000, 200, 0, 10);
        plots.AddHisto2F(selLabels_[i] + "_Ih_oldCorr_vs_RunNumber", 5000, 378000, 388000, 200, 0, 10);
        plots.AddHisto2F(selLabels_[i] + "_NPV_vs_RunNumber", 5000, 378000, 388000, 150, 0, 150);
        plots.AddHisto2F(selLabels_[i] + "_10000oP_vs_Ih", 300, 0, 300, 200, 0, 10);
        
        plots.AddHisto2F(selLabels_[i] + "_Fpix_vs_IhnearC", 200, C-0.3, C+0.3, 20, 0, 1);
        plots.AddHisto2F(selLabels_[i] + "_pT_vs_Fpixel", 50, 0, 1, 100, 0, 1000);

        plots.AddHisto1F(selLabels_[i] + "_electron_pt", 2000, 0, 2000);
        plots.AddHisto1F(selLabels_[i] + "_electron_eta", 60, -3, +3);
        plots.AddHisto1F(selLabels_[i] + "_electron_phi", 64, -3.2, 3.2);
        plots.AddHisto1F(selLabels_[i] + "_muon_pt", 2000, 0, 2000);
        plots.AddHisto1F(selLabels_[i] + "_muon_eta", 60, -3, +3);
        plots.AddHisto1F(selLabels_[i] + "_muon_phi", 64, -3.2, 3.2);

        // Syst on signal
        plots.AddHisto1F(selLabels_[i] + "_9fp10_SignalMass_nominal", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_9fp10_SignalMass_PUUp", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_9fp10_SignalMass_PUDown", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_9fp10_SignalMass_TriggerSFUp", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_9fp10_SignalMass_TriggerSFDown", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_9fp10_SignalMass_KUp", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_9fp10_SignalMass_KDown", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_9fp10_SignalMass_CUp", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_9fp10_SignalMass_CDown", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_9fp10_SignalMass_JetUp", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_9fp10_SignalMass_JetDown", 400, 0, 4000);

        plots.AddHisto1F(selLabels_[i] + "_9fp10_SignalMass_ATLASdata", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_9fp10_SignalMass_ATLASbkg", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_9fp10_SignalMass_ATLASglupion", 400, 0, 4000);

        plots.AddHisto1F(selLabels_[i] + "_99fp10_SignalMass_nominal", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_99fp10_SignalMass_PUUp", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_99fp10_SignalMass_PUDown", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_99fp10_SignalMass_TriggerSFUp", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_99fp10_SignalMass_TriggerSFDown", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_99fp10_SignalMass_KUp", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_99fp10_SignalMass_KDown", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_99fp10_SignalMass_CUp", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_99fp10_SignalMass_CDown", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_99fp10_SignalMass_JetUp", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_99fp10_SignalMass_JetDown", 400, 0, 4000);

        plots.AddHisto1F(selLabels_[i] + "_99fp10_SignalMass_ATLASbkg", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_99fp10_SignalMass_ATLASglupion", 400, 0, 4000);


        plots.AddHisto1F(selLabels_[i] + "_8fp9_SignalMass_nominal", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_8fp9_SignalMass_PUUp", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_8fp9_SignalMass_PUDown", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_8fp9_SignalMass_TriggerSFUp", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_8fp9_SignalMass_TriggerSFDown", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_8fp9_SignalMass_KUp", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_8fp9_SignalMass_KDown", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_8fp9_SignalMass_CUp", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_8fp9_SignalMass_CDown", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_8fp9_SignalMass_JetUp", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_8fp9_SignalMass_JetDown", 400, 0, 4000);

        plots.AddHisto1F(selLabels_[i] + "_8fp9_SignalMass_ATLASbkg", 400, 0, 4000);
        plots.AddHisto1F(selLabels_[i] + "_8fp9_SignalMass_ATLASglupion", 400, 0, 4000);



        vcp.push_back(std::move(plots));
    }

        // Without selections
    CPlots plots;
    plots.AddHisto1D("CandidateCutflow", 20, 0, 20);
    plots.AddHisto1D("EventCutflow", 20, 0, 20);
    plots.AddHisto1D("EventCutflow_NotrackCut", 20, 0, 20);
    plots.AddHisto1D("EventCutflow__isGen_lastbin", 2, -0.5, 1.5);

    plots.AddHisto1F("nHSCP", 10, 0, 10);
    plots.AddHisto1F("nPVgood", 150, 0, 150);
    plots.AddHisto1F("ndEdx_StripOnly", 50, 0, 50);

    plots.AddHisto1F("Trigger", 3, -1, 2);

    plots.AddHisto1F("Nm1_trigger", 4, -1.5, 2.5);
    plots.AddHisto1F("Nm1_METfilters", 4, -1.5, 2.5);
    plots.AddHisto1F("Nm1_PUppiMET", 200, 0, 2000);
    plots.AddHisto1F("Nm1_Ptpseudo", 300, 0, 3000);
    plots.AddHisto1F("Nm1_eta", 60, -3, +3);
    plots.AddHisto1F("Nm1_NOPH", 15, 0, 15);
    plots.AddHisto1F("Nm1_FOVH", 50, 0, 1.1);
    plots.AddHisto1F("Nm1_NOM", 50, 0, 50);
    plots.AddHisto1F("Nm1_HighPurity", 4, -1.5, 2.5);
    plots.AddHisto1F("Nm1_Chi2", 14, 0, 14);
    plots.AddHisto1F("Nm1_dZ", 100, -0.8, 0.8);
    plots.AddHisto1F("Nm1_dXY", 100, -0.1, 0.1);
    plots.AddHisto1F("Nm1_PFMiniIso", 1000, 0, 0.5);
    plots.AddHisto1F("Nm1_TrkIso", 400, 0, 200);
    plots.AddHisto1F("Nm1_EoverP", 500, 0, 5);
    plots.AddHisto1F("Nm1_PtErr_over_PtPt", 500, 0, 0.01);
    plots.AddHisto1F("Nm1_Fpix", 21, 0, 1.1);
    plots.AddHisto1F("Nm1_PtErr_over_Pt", 800, 0, 4);
    plots.AddHisto1F("Nm1_Ih_StripOnly", 200, 0, 10);
    plots.AddHisto1F("Nm1_Ih_StripOnly_rescaled", 200, 0, 10);


    plots.AddHisto1F("Nm1_event_trigger",           4,    -1.5,  2.5);
    plots.AddHisto1F("Nm1_event_METfilters",        4,    -1.5,  2.5);
    plots.AddHisto1F("Nm1_event_PUppiMET",           200,   0,    2000);
    plots.AddHisto1F("Nm1_event_Ptpseudo",          300,   0,    3000);
    plots.AddHisto1F("Nm1_event_eta",               60,   -3,    3);
    plots.AddHisto1F("Nm1_event_NOPH",              15,    0,    15);
    plots.AddHisto1F("Nm1_event_FOVH",              50,    0,    1.1);
    plots.AddHisto1F("Nm1_event_NOM",               50,    0,    50);
    plots.AddHisto1F("Nm1_event_HighPurity",        4,    -1.5,  2.5);
    plots.AddHisto1F("Nm1_event_Chi2",              14,    0,    14);
    plots.AddHisto1F("Nm1_event_dZ",                100,  -0.8,  0.8);
    plots.AddHisto1F("Nm1_event_dXY",               100,  -0.1,  0.1);
    plots.AddHisto1F("Nm1_event_PFMiniIso",         1000,  0,    0.5);
    plots.AddHisto1F("Nm1_event_TrkIso",            400,   0,    200);
    plots.AddHisto1F("Nm1_event_EoverP",            500,   0,    5);
    plots.AddHisto1F("Nm1_event_PtErr_over_PtPt",   500,   0,    0.01);
    plots.AddHisto1F("Nm1_event_Fpix",              21,    0,    1.1);
    plots.AddHisto1F("Nm1_event_PtErr_over_Pt",     800,   0,    4);
    plots.AddHisto1F("Nm1_event_Ih_StripOnly",      200,   0,    10);
    plots.AddHisto1F("Nm1_event_Ih_StripOnly_rescaled",      200,   0,    10);

    plots.AddHisto1F("RecoPFMET", 400, 0, 2000);


    plots.AddHisto2F("LastBinEventCutflow___GenPt_vs_trackPt", 100, 0, 2500, 100, 0, 2500);
    plots.AddHisto1F("LastBinEventCutflow___HSCP_type", 6, -0.5, 5.5);
    plots.AddHisto1F("LastBinEventCutflow___PF_type", 250, 0, 250);

    plots.AddHisto1F("PostTrigger_PseudoCaloMET_SF_rescaled", 200, 0, 2000);
    plots.AddHisto1F("PostTrigger_PseudoCaloMET_SF_NOTrescaled", 200, 0, 2000);

    plots.AddHisto1F("PostTrigger_PseudoCaloMET_rescaled", 200, 0, 2000);
    plots.AddHisto1F("PostTrigger_PseudoCaloMET_NOTrescaled", 200, 0, 2000);


    // General plots
    plots.AddHisto1F("Nosel_PseudoCaloMET", 200, 0, 2000);
    plots.AddHisto1F("Nosel_PseudoCaloMET_rescaled", 200, 0, 2000);
    plots.AddHisto1F("Nosel_Ptpseudo", 300, 0, 3000);
    plots.AddHisto1F("Nosel_P", 200, 0, 8000);
    plots.AddHisto1F("Noselbis_P", 200, 0, 8000);
    plots.AddHisto1F("Nosel_10000oP", 300, 0, 300);
    plots.AddHisto2F("Nosel_10000oP_vs_Eta", 300, 0, 300, 60, -3, +3);
    plots.AddHisto1F("Nosel_eta", 60, -3, +3);
    plots.AddHisto1F("Nosel_phi", 64, -3.2, 3.2);
    plots.AddHisto1F("Nosel_PFMiniIso", 1000, 0, 0.5);
    plots.AddHisto1F("Nosel_TrkIso", 400, 0, 200);
    plots.AddHisto1F("Nosel_EoverP", 500, 0, 5);
    plots.AddHisto1F("Nosel_PtErr_over_PtPt", 500, 0, 0.01);
    plots.AddHisto1F("Nosel_Fpix", 21, 0, 1.1);
    plots.AddHisto1F("Nosel_PtErr_over_Pt", 800, 0, 4);
    plots.AddHisto1F("Nosel_Ih", 600, 0, 30);
    plots.AddHisto1F("Nosel_Ih_rescaled", 200, 0, 10);
    plots.AddHisto1F("Nosel_Ih_oldCorr", 200, 0, 10);
    plots.AddHisto1F("Nosel_Ih_noSF", 200, 0, 10);
    plots.AddHisto1F("Nosel_PthatQCD", 2000, 0, 2000);

    plots.AddHisto1F("PostTrigger_P", 200, 0, 8000);
    plots.AddHisto1F("PostTriggerbis_P", 200, 0, 8000);
    plots.AddHisto1F("PostTrigger_PseudoCaloMETcut_P", 200, 0, 8000);
    plots.AddHisto1F("PostTrigger_PUppiMETcut_P", 200, 0, 8000);
    plots.AddHisto1F("PostTrigger_PseudoCaloMETcut_PUppiMETcut_P", 200, 0, 8000);
    plots.AddHisto1F("PostTriggerWsf_PseudoCaloMETcut_P", 200, 0, 8000);
    plots.AddHisto1F("PostTriggerWsf_PUppiMETcut_P", 200, 0, 8000);
    plots.AddHisto1F("PostTriggerWsf_PseudoCaloMETcut_PUppiMETcut_P", 200, 0, 8000);
    plots.AddHisto1F("PostTriggerWsfup_PseudoCaloMETcut_P", 200, 0, 8000);
    plots.AddHisto1F("PostTriggerWsfup_PUppiMETcut_P", 200, 0, 8000);
    plots.AddHisto1F("PostTriggerWsfup_PseudoCaloMETcut_PUppiMETcut_P", 200, 0, 8000);
    plots.AddHisto1F("PostTriggerWsfdown_PseudoCaloMETcut_P", 200, 0, 8000);
    plots.AddHisto1F("PostTriggerWsfdown_PUppiMETcut_P", 200, 0, 8000);
    plots.AddHisto1F("PostTriggerWsfdown_PseudoCaloMETcut_PUppiMETcut_P", 200, 0, 8000);

    plots.AddHisto1F("PUppiMET__postTrigger", 100, 0, 2500);
    plots.AddHisto1F("PUppiMET__postTrigger__SF", 100, 0, 2500);

    plots.AddHisto2F("Nosel_Fpix_vs_RunNumber", 5000, 378000, 388000, 20, 0, 1);
    plots.AddHisto2F("Nosel_Ih_vs_RunNumber", 5000, 378000, 388000, 200, 0, 10);
    plots.AddHisto2F("Nosel_Ih_oldCorr_vs_RunNumber", 5000, 378000, 388000, 200, 0, 10);
    plots.AddHisto2F("Nosel_NPV_vs_RunNumber", 5000, 378000, 388000, 150, 0, 150);
    plots.AddHisto2F("Nosel_10000oP_vs_Ih", 300, 0, 300, 200, 0, 10);

    plots.AddHisto2F("Nosel_Fpix_vs_IhnearC", 200, C-0.3, C+0.3, 20, 0, 1);
    plots.AddHisto2F("Nosel_pT_vs_Fpixel", 50, 0, 1, 100, 0, 1000);
    plots.AddHisto2F("Nosel_PseudoMET_vs_PFMET", 200, 0, 2000, 200, 0, 2000);
    plots.AddHisto2F("Nosel_PUppiMETNoMu_VS_PseudoMET", 100, 0, 2500, 100, 0, 2500);
    plots.AddHisto2F("HSCPPartialsel_PUppiMET_VS_PseudoMET", 100, 0, 2500, 100, 0, 2500);
    plots.AddHisto2F("HSCPPartialsel_PUppiMETNoMu_VS_PseudoMET", 100, 0, 2500, 100, 0, 2500);


    plots.AddHisto2F("Nosel_PseudoMET_JetSize", 100, 0, 2500, 200, 0, 200);
    plots.AddHisto2F("Nosel_PUppiMET_JetSize", 100, 0, 2500, 200, 0, 200);    

    plots.AddHisto1F("Nosel_PUppiMET", 100, 0, 2500);
    plots.AddHisto1F("Nosel_PUppiMETup", 100, 0, 2500);
    plots.AddHisto1F("Nosel_PUppiMETdown", 100, 0, 2500);


    plots.AddHisto2F("trackPT_vs_trackPseudoTrackPT__PFmuon", 100, 0, 2500, 100, 0, 2500);
    plots.AddHisto2F("trackPT_vs_trackPseudoTrackPT__PFpion", 100, 0, 2500, 100, 0, 2500);
    plots.AddHisto2F("trackPT_vs_trackPseudoTrackPT__HSCPmatched", 100, 0, 5000, 100, 0, 5000);
    plots.AddHisto2F("trackPT_vs_trackPseudoTrackPT", 100, 0, 2500, 100, 0, 2500);
    plots.AddHisto2F("trackETA_vs_trackPseudoTrackETA", 60, -3, +3, 60, -3, +3);
    plots.AddHisto2F("trackPHI_vs_trackPseudoTrackPHI", 64, -3.2, 3.2, 64, -3.2, 3.2);

    plots.AddHisto1F("Nosel_GenHSCPcharged_Eta", 60, -3, +3);
    plots.AddHisto1F("HSCPsel_GenHSCPcharged_Eta", 60, -3, +3);

    plots.AddHisto1F("GenHSCPmatching__PFType", 500, 0, 500);
    plots.AddHisto1F("GenHSCPmatching__HSCP_type", 6, -0.5, 5.5);
    plots.AddHisto1F("GenHSCPmatching__PFType__if_ORtrigger", 500, 0, 500);

    plots.AddHisto1F("Nosel_Gen__PairType", 5, 0, 5);
    plots.AddHisto1F("Nosel_Gen__nHSCP", 10, 0, 10);
    plots.AddHisto1F("Nosel_Gen__nMatchedTrack", 10, 0, 10);


    plots.AddHisto1F("Nosel_GenHSCP_beta", 40, 0, 1);
    plots.AddHisto1F("Nosel_GenHSCP_beta_ifHLTMu50", 40, 0, 1);
    
    

    plots.AddHisto1F("Nosel_Gen__PairType__if_HLT_PFMET120_PFMHT120_IDTight", 5, 0, 5);
    plots.AddHisto1F("Nosel_Gen__PairType__if_HLT_PFHT500_PFMET100_PFMHT100_IDTight", 5, 0, 5);
    plots.AddHisto1F("Nosel_Gen__PairType__if_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", 5, 0, 5);
    plots.AddHisto1F("Nosel_Gen__PairType__if_HLT_MET105_IsoTrk50", 5, 0, 5);
    plots.AddHisto1F("Nosel_Gen__PairType__if_ORtrigger", 5, 0, 5);
    plots.AddHisto1F("Nosel_Gen__PairType__if_HLT_Mu50", 5, 0, 5);
    plots.AddHisto1F("Nosel_GenHSCPmatching__PairType__1HSCP_if_HLT_PFMET120_PFMHT120_IDTight", 5, 0, 5);
    plots.AddHisto1F("Nosel_GenHSCPmatching__PairType__1HSCP_if_HLT_PFHT500_PFMET100_PFMHT100_IDTight", 5, 0, 5);
    plots.AddHisto1F("Nosel_GenHSCPmatching__PairType__1HSCP_if_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", 5, 0, 5);
    plots.AddHisto1F("Nosel_GenHSCPmatching__PairType__1HSCP_if_HLT_MET105_IsoTrk50", 5, 0, 5);
    plots.AddHisto1F("Nosel_GenHSCPmatching__PairType__1HSCP_if_ORtrigger", 5, 0, 5);
    plots.AddHisto1F("Nosel_GenHSCPmatching__PairType__2HSCP_if_HLT_PFMET120_PFMHT120_IDTight", 5, 0, 5);
    plots.AddHisto1F("Nosel_GenHSCPmatching__PairType__2HSCP_if_HLT_PFHT500_PFMET100_PFMHT100_IDTight", 5, 0, 5);
    plots.AddHisto1F("Nosel_GenHSCPmatching__PairType__2HSCP_if_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", 5, 0, 5);
    plots.AddHisto1F("Nosel_GenHSCPmatching__PairType__2HSCP_if_HLT_MET105_IsoTrk50", 5, 0, 5);
    plots.AddHisto1F("Nosel_GenHSCPmatching__PairType__2HSCP_if_ORtrigger", 5, 0, 5);
    plots.AddHisto1F("Nosel_GenHSCPmatching__PairType__0HSCP_if_HLT_PFMET120_PFMHT120_IDTight", 5, 0, 5);
    plots.AddHisto1F("Nosel_GenHSCPmatching__PairType__0HSCP_if_HLT_PFHT500_PFMET100_PFMHT100_IDTight", 5, 0, 5);
    plots.AddHisto1F("Nosel_GenHSCPmatching__PairType__0HSCP_if_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", 5, 0, 5);
    plots.AddHisto1F("Nosel_GenHSCPmatching__PairType__0HSCP_if_HLT_MET105_IsoTrk50", 5, 0, 5);
    plots.AddHisto1F("Nosel_GenHSCPmatching__PairType__0HSCP_if_ORtrigger", 5, 0, 5);

    plots.AddHisto1F("Nosel_GenHSCPmatching__PairType__0HSCP", 5, 0, 5);
    plots.AddHisto1F("Nosel_GenHSCPmatching__PairType__1HSCP", 5, 0, 5);
    plots.AddHisto1F("Nosel_GenHSCPmatching__PairType__2HSCP", 5, 0, 5);



    plots.AddHisto1F("Nosel_PU_up", 100, 0, 100);
    plots.AddHisto1F("Nosel_PU_down", 100, 0, 100);
    plots.AddHisto1F("Nosel_PU", 100, 0, 100);
    plots.AddHisto1F("Nosel_PU_Factor_up", 100, 0, 10);
    plots.AddHisto1F("Nosel_PU_Factor_down", 100, 0, 10);


    plots.AddHisto2F("dEdX0stripVsP", 100, 0, 4000, 50, 0., 50.);
    plots.AddHisto2F("dEdX0stripVsP__large", 1000, 0, 4000, 500, 0.,50.);

    plots.AddHisto2F("Nosel_PUppiMET_VS_PseudoMET", 100, 0, 2500, 100, 0, 2500);

    
    
    vcp_nosel.push_back(std::move(plots));

    std::cout << std::endl;
    std::cout << "Histograms initialised" <<std::endl;
}

Bool_t HSCPSelector::Process(Long64_t entry) {
    fReader.SetLocalEntry(entry);

    //Step 1: Calculate PFHT, PFMHT and PFMHTNoMu variables for the event
    //        with associated up/nominal/down variations
    double PFHT = 0., PFHT_up = 0, PFHT_down = 0;
    double PFMHT_x = 0., PFMHT_x_up = 0, PFMHT_x_down = 0;
    double PFMHT_y = 0., PFMHT_y_up = 0, PFMHT_y_down = 0;

    double PUppiMET_x = RecoPUppiMET[0] * cos(RecoPUppiMET_phi[0]);
    double PUppiMET_y = RecoPUppiMET[0] * sin(RecoPUppiMET_phi[0]);
    double PUppiMET_x_up = PUppiMET_x,   PUppiMET_y_up = PUppiMET_y;
    double PUppiMET_x_down = PUppiMET_x, PUppiMET_y_down = PUppiMET_y;

    for ( unsigned int j = 0; j < Jet_pt.GetSize(); j++ ) {
        if (Jet_passJetID[j] && !Jet_IsRejectedbyJVM[j]) {
            double d = Jet_jesUncTotal[j];

            // nominal
            PFHT    += std::abs(Jet_pt[j]);
            PFMHT_x += Jet_px[j];
            PFMHT_y += Jet_py[j];

            // JES up
            PFHT_up    += std::abs(Jet_pt[j]) * (1. + d);
            PFMHT_x_up += Jet_px[j] * (1. + d);
            PFMHT_y_up += Jet_py[j] * (1. + d);

            // JES down
            PFHT_down    += std::abs(Jet_pt[j]) * (1. - d);
            PFMHT_x_down += Jet_px[j] * (1. - d);
            PFMHT_y_down += Jet_py[j] * (1. - d);

            // PUppiMET up/down
            PUppiMET_x_up   -= Jet_px[j] * d;
            PUppiMET_y_up   -= Jet_py[j] * d;

            PUppiMET_x_down += Jet_px[j] * d;
            PUppiMET_y_down += Jet_py[j] * d;
        }
    }

    double PUppiMET_up   = sqrt(PUppiMET_x_up*PUppiMET_x_up + PUppiMET_y_up*PUppiMET_y_up);
    double PUppiMET_down = sqrt(PUppiMET_x_down*PUppiMET_x_down + PUppiMET_y_down*PUppiMET_y_down);

    double PFMHTNoMu_x = PFMHT_x, PFMHTNoMu_x_up = PFMHT_x_up, PFMHTNoMu_x_down = PFMHT_x_down;
    double PFMHTNoMu_y = PFMHT_y, PFMHTNoMu_y_up = PFMHT_y_up, PFMHTNoMu_y_down = PFMHT_y_down;

    for (unsigned m = 0; m <muon_pt.GetSize(); m++) {
        if (muon_isPFMuon[m]) {
            PFMHTNoMu_x -= muon_pt[m] * cos(muon_phi[m]);
            PFMHTNoMu_y -= muon_pt[m] * sin(muon_phi[m]);

            PFMHTNoMu_x_up -= muon_pt[m] * cos(muon_phi[m]);
            PFMHTNoMu_y_up -= muon_pt[m] * sin(muon_phi[m]);

            PFMHTNoMu_x_down -= muon_pt[m] * cos(muon_phi[m]);
            PFMHTNoMu_y_down -= muon_pt[m] * sin(muon_phi[m]);

            // PUppiMET NoMu
            PUppiMET_x -= muon_pt[m] * cos(muon_phi[m]);
            PUppiMET_y -= muon_pt[m] * sin(muon_phi[m]);

            PUppiMET_x_up -= muon_pt[m] * cos(muon_phi[m]);
            PUppiMET_y_up -= muon_pt[m] * sin(muon_phi[m]);
            
            PUppiMET_x_down -= muon_pt[m] * cos(muon_phi[m]);
            PUppiMET_y_down -= muon_pt[m] * sin(muon_phi[m]);
        }
    }

    double PFMHT      = sqrt(PFMHT_x*PFMHT_x + PFMHT_y*PFMHT_y);
    double PFMHT_up   = sqrt(PFMHT_x_up*PFMHT_x_up + PFMHT_y_up*PFMHT_y_up);
    double PFMHT_down = sqrt(PFMHT_x_down*PFMHT_x_down + PFMHT_y_down*PFMHT_y_down);

    double PFMHTNoMu      = sqrt(PFMHTNoMu_x*PFMHTNoMu_x + PFMHTNoMu_y*PFMHTNoMu_y);
    double PFMHTNoMu_up   = sqrt(PFMHTNoMu_x_up*PFMHTNoMu_x_up + PFMHTNoMu_y_up*PFMHTNoMu_y_up);
    double PFMHTNoMu_down = sqrt(PFMHTNoMu_x_down*PFMHTNoMu_x_down + PFMHTNoMu_y_down*PFMHTNoMu_y_down);
    
    double PUppiMET_NoMu      = sqrt(PUppiMET_x*PUppiMET_x + PUppiMET_y*PUppiMET_y);

    
    //----------------------------------
    //Loop over all HSCP candidates
    //----------------------------------
    vector<int> iCand(selLabels_.size(),-1);
    vector<float> maxIh(selLabels_.size(),-1);

    std::vector<bool> eventCuts(19, false);
    std::vector<bool> eventCuts_noTrackCut(19, false);
    std::vector<bool> singleCut(19, false);
    std::vector<bool> eventPassedNm1(19, false);
    std::vector<unsigned int> iCandNm1(19, (int)-1);
    std::vector<float> maxIhNm1(19, -999.);

    bool trigger = false, triggerONLY = false, trigger_varUP = false, trigger_varDOWN = false;
    bool METfilters = false;

    bool isMC = false;
    bool isSignal = (sample_.find("Gluino") != std::string::npos) || 
                    (sample_.find("Stop") != std::string::npos) || 
                    (sample_.find("Stau") != std::string::npos) ? true : false;

    float AppliedWeight = 1.;
    float AppliedWeight_SF_Up = 1., AppliedWeight_SF_Down = 1.;
    float AppliedWeight_PU_Up = 1., AppliedWeight_PU_Down = 1.;

    float AppliedWeight1Dpseudo = 1., AppliedWeight1Dpseudo_SF_Down = 1., AppliedWeight1Dpseudo_SF_Up = 1.;
    float AppliedWeight1Dpuppi = 1., AppliedWeight1Dpuppi_SF_Down = 1., AppliedWeight1Dpuppi_SF_Up = 1.;
    if (sample_.find("Mu2024") != std::string::npos ||
        sample_.find("TestMuon2024") != std::string::npos) {
        trigger = *HLT_Mu50;
        METfilters = true; // no MET filters in muon datasets
    }
    else if (sample_.find("MET") != std::string::npos || 
             sample_.find("QCD2024") != std::string::npos || 
             sample_.find("TTbar2024") != std::string::npos ||
             sample_.find("TTbarSemiLep2024") != std::string::npos ||
             sample_.find("Wjets") != std::string::npos ||
             sample_.find("Gluino") != std::string::npos ||
             sample_.find("Stop") != std::string::npos || 
             sample_.find("Stau") != std::string::npos ||
             sample_.find("WjetMuNu") != std::string::npos) {
        
        triggerONLY     = ( *HLT_PFMET120_PFMHT120_IDTight || *HLT_PFHT500_PFMET100_PFMHT100_IDTight ||
                        *HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 || *HLT_MET105_IsoTrk50);
        trigger         = ((*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) || (*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100)
                        || (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60) || *HLT_MET105_IsoTrk50);
        trigger_varUP   = ((*HLT_PFMET120_PFMHT120_IDTight && PFMHT_up>120) || (*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT_up>500 && PFMHT_up>100)
                        || (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu_up>120 && PFHT_up>60) || *HLT_MET105_IsoTrk50);
        trigger_varDOWN = ((*HLT_PFMET120_PFMHT120_IDTight && PFMHT_down>120) || (*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT_down>500 && PFMHT_down>100)
                        || (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu_down>120 && PFHT_down>60) || *HLT_MET105_IsoTrk50);

        METfilters = Flag_allMETFilters[0];
    }
    if (sample_.find("Gluino") != std::string::npos ||
        sample_.find("Stop") != std::string::npos || 
        sample_.find("Stau") != std::string::npos ||
        sample_.find("QCD2024") != std::string::npos || 
        sample_.find("TTbar2024") != std::string::npos ||
        sample_.find("TTbarSemiLep2024") != std::string::npos ||
        sample_.find("Wjets") != std::string::npos ||
        sample_.find("WjetMuNu") != std::string::npos) {
        isMC = true;


        float SF_Down = 1., SF = 1., SF_Up = 1.;
        bool found = false;
        for (size_t i = 0; i < SF2D_PseudoMETlovalue.size(); ++i) {
            if (PseudoCaloMET[0] >= SF2D_PseudoMETlovalue[i] && PseudoCaloMET[0] < SF2D_PseudoMEThivalue[i] &&
                RecoPUppiMET[0]  >= SF2D_PUppiMETlovalue[i]  && RecoPUppiMET[0]  < SF2D_PUppiMEThivalue[i]) {
                SF_Down = SF2D_Down[i];
                SF      = SF2D[i];
                SF_Up   = SF2D_Up[i];
                found = true;
                break;
            }
        }
        if (!found) {
            // candidat hors couverture (ex. au-delà du dernier high edge, ou sentinel -2)
            SF_Down = SF2D_Down.back();
            SF      = SF2D.back();
            SF_Up   = SF2D_Up.back();
        }


        // 1D SF for PseudoMET
        float SF1Dpseudo_Down_temp = 1., SF1Dpseudo_temp = 1., SF1Dpseudo_Up_temp = 1.;
        for (size_t i = 0; i < SF1D_PseudoMETvalue.size(); ++i) {
            float hi = (i+1 < SF1D_PseudoMETvalue.size())
                        ? SF1D_PseudoMETvalue[i+1]
                        : std::numeric_limits<float>::max();
            if (PseudoCaloMET[0] >= SF1D_PseudoMETvalue[i] && PseudoCaloMET[0] < hi) {
                SF1Dpseudo_Down_temp = SF1Dpseudo_Down[i];
                SF1Dpseudo_temp      = SF1Dpseudo[i];
                SF1Dpseudo_Up_temp   = SF1Dpseudo_Up[i];
                break;
            }
        }

        // 1D SF for PUppiMET
        float SF1Dpuppi_Down_temp = 1., SF1Dpuppi_temp = 1., SF1Dpuppi_Up_temp = 1.;
        for (size_t i = 0; i < SF1D_PUppiMETvalue.size(); ++i) {
            float hi = (i+1 < SF1D_PUppiMETvalue.size())
                        ? SF1D_PUppiMETvalue[i+1]
                        : std::numeric_limits<float>::max();
            if (RecoPUppiMET[0] >= SF1D_PUppiMETvalue[i] && RecoPUppiMET[0] < hi) {
                SF1Dpuppi_Down_temp = SF1Dpuppi_Down[i];
                SF1Dpuppi_temp      = SF1Dpuppi[i];
                SF1Dpuppi_Up_temp   = SF1Dpuppi_Up[i];
                break;
            }
        }


        // 2D case
        AppliedWeight_SF_Down = *weightPU * SF_Down;
        AppliedWeight = *weightPU * SF;
        AppliedWeight_SF_Up = *weightPU * SF_Up;

        // 1D PseudoMET
        AppliedWeight1Dpseudo_SF_Down = *weightPU * SF1Dpseudo_Down_temp;
        AppliedWeight1Dpseudo = *weightPU * SF1Dpseudo_temp;
        AppliedWeight1Dpseudo_SF_Up = *weightPU * SF1Dpseudo_Up_temp;

        // 1D PUppiMET
        AppliedWeight1Dpuppi_SF_Down = *weightPU * SF1Dpuppi_Down_temp;
        AppliedWeight1Dpuppi = *weightPU * SF1Dpuppi_temp;
        AppliedWeight1Dpuppi_SF_Up = *weightPU * SF1Dpuppi_Up_temp;

        
        // PU up/down
        AppliedWeight_PU_Up = *weightPU_Up * *weightPU * SF1Dpuppi_temp;
        AppliedWeight_PU_Down = *weightPU_Down * *weightPU * SF1Dpuppi_temp;
    }



    singleCut[0] = trigger;
    if (singleCut[0]) {
        vcp_nosel[0].FillHisto1D("EventCutflow", 1.5, AppliedWeight1Dpuppi);
        vcp_nosel[0].FillHisto1D("EventCutflow_NotrackCut", 1.5, AppliedWeight1Dpuppi);
    }
    singleCut[1] = singleCut[0] && METfilters;
    if (singleCut[1]) {
        vcp_nosel[0].FillHisto1D("EventCutflow", 2.5, AppliedWeight1Dpuppi);
        vcp_nosel[0].FillHisto1D("EventCutflow_NotrackCut", 2.5, AppliedWeight1Dpuppi);
    }
    singleCut[2] = singleCut[1] && RecoPUppiMET[0] > 150.;
    if (singleCut[2]) {
        vcp_nosel[0].FillHisto1D("EventCutflow", 3.5, AppliedWeight1Dpuppi);
        vcp_nosel[0].FillHisto1D("EventCutflow_NotrackCut", 3.5, AppliedWeight1Dpuppi);
    }


    unsigned int i = 0;
    for (unsigned int j=0; j<HSCP_hasTrack.GetSize(); j++){

        if (!HSCP_hasTrack[j]) continue;

            //TAKE MOST IONIZING CANDIDATE
        for(unsigned int s=0;s<selections_.size();s++){
            bool (HSCPSelector::*ptr)(int);
            ptr = selections_[s];
            if((this->*ptr)(i) && (trigger)){// || trigger_varUP || trigger_varDOWN)){
                if(Ih_Strip[i] > maxIh[s]){
                    maxIh[s] = Ih_Strip[i];
                    iCand[s] = i;
                }
            }
        }


        // CUTFLOW
        singleCut[3] = singleCut[2] && (Pt_pseudo[i] > 50.) && (Pt[i] > 50.);
        singleCut[4] = singleCut[3] && fabs(Eta[i]) < 2.4;
        singleCut[5] = singleCut[4] && (NbPixelHit_noL1[i] >= 2);
        singleCut[6] = singleCut[5] && (NOM_noL1[i] >= 10);
        singleCut[7] = singleCut[6] && (FracOfValidHit[i] > 0.8);
        singleCut[8] = singleCut[7] && isHighPurityTrack[i];
        singleCut[9] = singleCut[8] && (normChi2[i] < 5);
        singleCut[10] = singleCut[9] && (fabs(dz[i]) < 0.1);
        singleCut[11] = singleCut[10] && (fabs(dxy[i]) < 0.02);
        singleCut[12] = singleCut[11] && (miniRelIsoAll[i] < 0.02);
        singleCut[13] = singleCut[12] && (IsoSumPt_dr03[i] < 15);
        singleCut[14] = singleCut[13] && (EoP[i] < 0.3);
        singleCut[15] = singleCut[14] && (ptOverptErrptErr[i] < 0.0008);
        singleCut[16] = singleCut[15] && (ptOverptErr[i] < 1);
        singleCut[17] = singleCut[16] && (Fpix[i] > 0.3);
        singleCut[18] = singleCut[17] && (Ih_Strip[i] > C);

            // candidate CUTFLOW
        vcp_nosel[0].FillHisto1D("CandidateCutflow", 0, *weightPU); // All candidates with a track
        if (singleCut[0]) vcp_nosel[0].FillHisto1D("CandidateCutflow", 1, AppliedWeight1Dpuppi);
        if (singleCut[1]) vcp_nosel[0].FillHisto1D("CandidateCutflow", 2, AppliedWeight1Dpuppi);
        for (unsigned int j = 2; j < singleCut.size(); j++) {
            if (singleCut[j]) vcp_nosel[0].FillHisto1D("CandidateCutflow", j+1, AppliedWeight1Dpuppi);
        }


        std::vector<std::function<bool(int)>> cuts;
        cuts.push_back([&](int i){ return trigger; });
        cuts.push_back([&](int i){ return METfilters; });
        cuts.push_back([&](int i){ return RecoPUppiMET[0] > 150.; });
        cuts.push_back([&](int i){ return ((Pt_pseudo[i] > 50.0) && (Pt[i] > 50.0)); });
        cuts.push_back([&](int i){ return fabs(Eta[i]) < 2.4; });
        cuts.push_back([&](int i){ return NbPixelHit_noL1[i] >= 2; });
        cuts.push_back([&](int i){ return NOM_noL1[i] >= 10; });
        cuts.push_back([&](int i){ return FracOfValidHit[i] > 0.8; });
        cuts.push_back([&](int i){ return isHighPurityTrack[i]; });
        cuts.push_back([&](int i){ return normChi2[i] < 5.0; });
        cuts.push_back([&](int i){ return fabs(dz[i]) < 0.1; });
        cuts.push_back([&](int i){ return fabs(dxy[i]) < 0.02; });
        cuts.push_back([&](int i){ return miniRelIsoAll[i] < 0.02; });
        cuts.push_back([&](int i){ return IsoSumPt_dr03[i] < 15; });
        cuts.push_back([&](int i){ return EoP[i] < 0.3; });
        cuts.push_back([&](int i){ return ptOverptErrptErr[i] < 0.0008; });
        cuts.push_back([&](int i){ return ptOverptErr[i] < 1; });
        cuts.push_back([&](int i){ return Fpix[i] > 0.3; });
        cuts.push_back([&](int i){ return Ih_Strip[i] > C; });
        

            // Event CUTFLOW setup
        bool allPassed = true;
        for (unsigned int j = 0; j < cuts.size(); ++j) {
            if (allPassed && cuts[j](i)) eventCuts[j] = true;
            else allPassed = false;
        }

            // Event CUTFLOW without track cuts
        std::vector<std::function<bool(int)>> cuts_NotrackCut = cuts;
        cuts_NotrackCut[3] = [&](int i){ return true; }; // Pt
        cuts_NotrackCut[12] = [&](int i){ return true; }; // miniRelIsoAll
        cuts_NotrackCut[13] = [&](int i){ return true; }; // IsoSumPt_dr03
        cuts_NotrackCut[14] = [&](int i){ return true; }; // EoP
        cuts_NotrackCut[15] = [&](int i){ return true; }; // ptOverptErrptErr
        cuts_NotrackCut[17] = [&](int i){ return true; }; // ptOverptErr

        bool allPassed_NotrackCut = true;
        for (unsigned int j = 0; j < cuts.size(); ++j) {
            if (allPassed_NotrackCut && cuts_NotrackCut[j](i)) eventCuts_noTrackCut[j] = true;
            else allPassed_NotrackCut = false;
        }

            // N-1 distributions
        std::vector <bool> passedCuts(cuts.size(), true);
        for (unsigned int icut = 0; icut < cuts.size(); ++icut) {
            bool passedAllButOne = true;
            for (unsigned int jcut = 0; jcut < cuts.size(); ++jcut) {
                if (jcut == icut) continue;
                if (!cuts[jcut](i)) { passedAllButOne = false; break; }
            }
            passedCuts[icut] = passedAllButOne; // N-1 candidat
            if (passedAllButOne) {              // N-1 event
                eventPassedNm1[icut] = true;
                if (Ih_Strip[i] > maxIhNm1[icut]) {
                    maxIhNm1[icut] = Ih_Strip[i];
                    iCandNm1[icut] = i;
                }
            }
        }


        if (passedCuts[0]) vcp_nosel[0].FillHisto1F("Nm1_trigger", trigger, AppliedWeight1Dpuppi);
        if (passedCuts[1]) vcp_nosel[0].FillHisto1F("Nm1_METfilters", METfilters, AppliedWeight1Dpuppi);
        if (passedCuts[2]) vcp_nosel[0].FillHisto1F("Nm1_PUppiMET", RecoPUppiMET[0], AppliedWeight1Dpuppi);
        if (passedCuts[3]) vcp_nosel[0].FillHisto1F("Nm1_Ptpseudo", Pt_pseudo[i], AppliedWeight1Dpuppi);
        if (passedCuts[4]) vcp_nosel[0].FillHisto1F("Nm1_eta", Eta[i], AppliedWeight1Dpuppi);
        if (passedCuts[5]) vcp_nosel[0].FillHisto1F("Nm1_NOPH", NbPixelHit_noL1[i], AppliedWeight1Dpuppi);
        if (passedCuts[6]) vcp_nosel[0].FillHisto1F("Nm1_NOM", NOM_noL1[i], AppliedWeight1Dpuppi);
        if (passedCuts[7]) vcp_nosel[0].FillHisto1F("Nm1_FOVH", FracOfValidHit[i], AppliedWeight1Dpuppi);
        if (passedCuts[8]) vcp_nosel[0].FillHisto1F("Nm1_HighPurity", isHighPurityTrack[i], AppliedWeight1Dpuppi);
        if (passedCuts[9]) vcp_nosel[0].FillHisto1F("Nm1_Chi2", normChi2[i], AppliedWeight1Dpuppi);
        if (passedCuts[10]) vcp_nosel[0].FillHisto1F("Nm1_dZ", dz[i], AppliedWeight1Dpuppi);
        if (passedCuts[11]) vcp_nosel[0].FillHisto1F("Nm1_dXY", dxy[i], AppliedWeight1Dpuppi);
        if (passedCuts[12]) vcp_nosel[0].FillHisto1F("Nm1_PFMiniIso", miniRelIsoAll[i], AppliedWeight1Dpuppi);
        if (passedCuts[13]) vcp_nosel[0].FillHisto1F("Nm1_TrkIso", IsoSumPt_dr03[i], AppliedWeight1Dpuppi);
        if (passedCuts[14]) vcp_nosel[0].FillHisto1F("Nm1_EoverP", EoP[i], AppliedWeight1Dpuppi);
        if (passedCuts[15]) vcp_nosel[0].FillHisto1F("Nm1_PtErr_over_PtPt", ptOverptErrptErr[i], AppliedWeight1Dpuppi);
        if (passedCuts[16]) vcp_nosel[0].FillHisto1F("Nm1_PtErr_over_Pt", ptOverptErr[i], AppliedWeight1Dpuppi);
        if (passedCuts[17] && isMC) vcp_nosel[0].FillHisto1F("Nm1_Fpix", Fpix[i], AppliedWeight1Dpuppi);
        else if (passedCuts[17] && !isMC && Fpix[i] <= 0.9) vcp_nosel[0].FillHisto1F("Nm1_Fpix", Fpix[i], AppliedWeight1Dpuppi);
        if (passedCuts[18] && isSignal) vcp_nosel[0].FillHisto1F("Nm1_Ih_StripOnly", Ih_Strip[i], AppliedWeight1Dpuppi);
        else if (passedCuts[18] && !isSignal && Fpix[i] <= 0.9) vcp_nosel[0].FillHisto1F("Nm1_Ih_StripOnly", Ih_Strip[i], AppliedWeight1Dpuppi);
        if (passedCuts[18] && isSignal) vcp_nosel[0].FillHisto1F("Nm1_Ih_StripOnly_rescaled", Ih_Strip[i]*3.2168/3.12692, AppliedWeight1Dpuppi);       // MPV_data/MPV_MC
        else if (passedCuts[18] && !isSignal && Fpix[i] <= 0.9) vcp_nosel[0].FillHisto1F("Nm1_Ih_StripOnly_rescaled", Ih_Strip[i]*3.2168/3.12692, AppliedWeight1Dpuppi);


            // Event cutflow last bin: check GenPart info and PF_type
        if (eventCuts[18]) {
            for (unsigned int j = 0; j < GenPart_pt.GetSize(); j++) {
                float dEta = GenPart_eta[j] - Eta[i];
                float dPhi = GenPart_phi[j] - Phi[i];
                while (dPhi >  M_PI) dPhi -= 2*M_PI;
                while (dPhi < -M_PI) dPhi += 2*M_PI;
                if (std::sqrt(dEta*dEta + dPhi*dPhi) < 0.01) {
                    vcp_nosel[0].FillHisto2F("LastBinEventCutflow___GenPt_vs_trackPt", GenPart_pt[j], Pt_pseudo[i], AppliedWeight1Dpuppi);
                    vcp_nosel[0].FillHisto1D("EventCutflow__isGen_lastbin", 1, AppliedWeight1Dpuppi);
                }
            }
            vcp_nosel[0].FillHisto1F("LastBinEventCutflow___HSCP_type", HSCP_type[i], AppliedWeight1Dpuppi);    // No matching on the track + no PDG id check ! Careful please
            vcp_nosel[0].FillHisto1F("LastBinEventCutflow___PF_type", PF_type[i], AppliedWeight1Dpuppi);
        }

        i++;
    } //End of loop over all HSCP candidates

    // Event CUTFLOW
    vcp_nosel[0].FillHisto1D("EventCutflow", 0.5, *weightPU); // All events
    for (unsigned int j = 3; j < 19; ++j) { // cuts length
        if (eventCuts[j]) vcp_nosel[0].FillHisto1D("EventCutflow", j+1.5, AppliedWeight1Dpuppi);
    }

    vcp_nosel[0].FillHisto1D("EventCutflow_NotrackCut", 0.5, AppliedWeight1Dpuppi); // All events
    for (unsigned int j = 3; j < 19; ++j) { // cuts length
        if (eventCuts_noTrackCut[j]) vcp_nosel[0].FillHisto1D("EventCutflow_NotrackCut", j+1.5, AppliedWeight1Dpuppi);
    }


    // N-1 EVENT CUTFLOW histograms
    if (eventPassedNm1[0])  vcp_nosel[0].FillHisto1F("Nm1_event_trigger",          trigger, AppliedWeight1Dpuppi);
    if (eventPassedNm1[1])  vcp_nosel[0].FillHisto1F("Nm1_event_METfilters",       METfilters, AppliedWeight1Dpuppi);
    if (eventPassedNm1[2])  vcp_nosel[0].FillHisto1F("Nm1_event_PUppiMET",         RecoPUppiMET[0], AppliedWeight1Dpuppi);
    if (eventPassedNm1[3])  vcp_nosel[0].FillHisto1F("Nm1_event_Ptpseudo",         Pt_pseudo[iCandNm1[3]], AppliedWeight1Dpuppi);
    if (eventPassedNm1[4])  vcp_nosel[0].FillHisto1F("Nm1_event_eta",              Eta[iCandNm1[4]], AppliedWeight1Dpuppi);
    if (eventPassedNm1[5])  vcp_nosel[0].FillHisto1F("Nm1_event_NOPH",             NbPixelHit_noL1[iCandNm1[5]], AppliedWeight1Dpuppi);
    if (eventPassedNm1[6])  vcp_nosel[0].FillHisto1F("Nm1_event_NOM",              NOM_noL1[iCandNm1[6]], AppliedWeight1Dpuppi);
    if (eventPassedNm1[7])  vcp_nosel[0].FillHisto1F("Nm1_event_FOVH",             FracOfValidHit[iCandNm1[7]], AppliedWeight1Dpuppi);
    if (eventPassedNm1[8])  vcp_nosel[0].FillHisto1F("Nm1_event_HighPurity",       isHighPurityTrack[iCandNm1[8]], AppliedWeight1Dpuppi);
    if (eventPassedNm1[9])  vcp_nosel[0].FillHisto1F("Nm1_event_Chi2",             normChi2[iCandNm1[9]], AppliedWeight1Dpuppi);
    if (eventPassedNm1[10]) vcp_nosel[0].FillHisto1F("Nm1_event_dZ",               dz[iCandNm1[10]], AppliedWeight1Dpuppi);
    if (eventPassedNm1[11]) vcp_nosel[0].FillHisto1F("Nm1_event_dXY",              dxy[iCandNm1[11]], AppliedWeight1Dpuppi);
    if (eventPassedNm1[12]) vcp_nosel[0].FillHisto1F("Nm1_event_PFMiniIso",        miniRelIsoAll[iCandNm1[12]], AppliedWeight1Dpuppi);
    if (eventPassedNm1[13]) vcp_nosel[0].FillHisto1F("Nm1_event_TrkIso",           IsoSumPt_dr03[iCandNm1[13]], AppliedWeight1Dpuppi);
    if (eventPassedNm1[14]) vcp_nosel[0].FillHisto1F("Nm1_event_EoverP",           EoP[iCandNm1[14]], AppliedWeight1Dpuppi);
    if (eventPassedNm1[15]) vcp_nosel[0].FillHisto1F("Nm1_event_PtErr_over_PtPt",  ptOverptErrptErr[iCandNm1[15]], AppliedWeight1Dpuppi);
    if (eventPassedNm1[16]) vcp_nosel[0].FillHisto1F("Nm1_event_PtErr_over_Pt",    ptOverptErr[iCandNm1[16]], AppliedWeight1Dpuppi);
    if (eventPassedNm1[17] && isMC)  vcp_nosel[0].FillHisto1F("Nm1_event_Fpix",    Fpix[iCandNm1[17]], AppliedWeight1Dpuppi);
    else if (eventPassedNm1[17] && !isMC && Fpix[iCandNm1[17]] <= 0.9) vcp_nosel[0].FillHisto1F("Nm1_event_Fpix", Fpix[iCandNm1[17]], AppliedWeight1Dpuppi);
    if (eventPassedNm1[18] && isSignal) vcp_nosel[0].FillHisto1F("Nm1_event_Ih_StripOnly", Ih_Strip[iCandNm1[18]], AppliedWeight1Dpuppi);
    else if (eventPassedNm1[18] && !isSignal && Fpix[iCandNm1[18]] <= 0.9) vcp_nosel[0].FillHisto1F("Nm1_event_Ih_StripOnly", Ih_Strip[iCandNm1[18]], AppliedWeight1Dpuppi);
    if (eventPassedNm1[18] && isSignal) vcp_nosel[0].FillHisto1F("Nm1_event_Ih_StripOnly_rescaled", Ih_Strip[iCandNm1[18]]*3.2168/3.12692, AppliedWeight1Dpuppi);
    else if (eventPassedNm1[18] && !isSignal && Fpix[iCandNm1[18]] <= 0.9) vcp_nosel[0].FillHisto1F("Nm1_event_Ih_StripOnly_rescaled", Ih_Strip[iCandNm1[18]]*3.2168/3.12692, AppliedWeight1Dpuppi);



    // PseudoCaloMET variable calibration
    for(unsigned int s=0;s<selections_.size();s++) {
        if (selLabels_[s] == "CalibPseudoMET" || selLabels_[s] == "CalibPseudoMET_isRescaled") {
            bool BasicSel = false;
            unsigned int hasPassedMuon = 0;
            unsigned int hasPassedElectron = 0;
            unsigned int hadPassedHSCP = 0;

            if (*HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ && 
                *HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL && 
                muon_pt.GetSize()>0 && electron_pt.GetSize()>0 && Flag_allMETFilters[0]) BasicSel = true;

            for (unsigned int im = 0; im < muon_pt.GetSize(); im++) {
                if (muon_isTight[im] && muon_pt[im]>50 && muon_pfMiniRelIsoAll[im]<0.15) hasPassedMuon++;
            }
            for (unsigned int ie = 0; ie < electron_pt.GetSize(); ie++) {
                if (electron_isTight[ie]==1 && electron_pt[ie]>20 && electron_pfMiniRelIsoAll[ie]<0.15) hasPassedElectron++;
            }

            unsigned int i_track = 0;
            for(unsigned int j=0; j<HSCP_hasTrack.GetSize(); j++){    

                if (!HSCP_hasTrack[j]) continue;

                if ( (Pt[i_track] > 50.0) && (Pt_pseudo[i_track] > 50.0) && (std::abs(Eta[i_track]) < 2.4) && (NbPixelHit_noL1[i_track] >= 2) && (FracOfValidHit[i_track] > 0.8) && 
                (NOM_noL1[i_track] >= 10) && (isHighPurityTrack[i_track] == true) && (normChi2[i_track] < 5.0) && (std::abs(dz[i_track]) < 0.1) && (std::abs(dxy[i_track]) < 0.02) && 
                (miniRelIsoAll[i_track] < 0.02) && (IsoSumPt_dr03[i_track] < 15) && (ptOverptErrptErr[i_track] < 0.0008) ) hadPassedHSCP++;

                i_track++;
            }
            
            if (BasicSel && hasPassedMuon==1 && hasPassedElectron==1 && hadPassedHSCP>0) {

                float isRescaled = (selLabels_[s] == "CalibPseudoMET_isRescaled") ? 163.451/146.864 : 1;
                float weightOnCalib = *weightPU;

                for (unsigned int im = 0; im < muon_pt.GetSize(); im++) {
                    if (muon_isTight[im] && muon_pt[im]>50 && muon_pfMiniRelIsoAll[im]<0.15) {
                        vcp[s].FillHisto1F(selLabels_[s] + "_muon_pt", muon_pt[im], weightOnCalib);
                        vcp[s].FillHisto1F(selLabels_[s] + "_muon_eta", muon_eta[im], weightOnCalib);
                        vcp[s].FillHisto1F(selLabels_[s] + "_muon_phi", muon_phi[im], weightOnCalib);
                    }
                }
                for (unsigned int ie = 0; ie < electron_pt.GetSize(); ie++) {
                    if (electron_isTight[ie]==1 && electron_pt[ie]>20 && electron_pfMiniRelIsoAll[ie]<0.15) {
                        vcp[s].FillHisto1F(selLabels_[s] + "_electron_pt", electron_pt[ie], weightOnCalib);
                        vcp[s].FillHisto1F(selLabels_[s] + "_electron_eta", electron_eta[ie], weightOnCalib);
                        vcp[s].FillHisto1F(selLabels_[s] + "_electron_phi", electron_phi[ie], weightOnCalib);
                    }
                }
                
                vcp[s].FillHisto1F(selLabels_[s] + "_PUppiMET", RecoPUppiMET[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_RecoPFMET", RecoPFMET[0], weightOnCalib);
                if (PseudoCaloMET[0] > 170) vcp[s].FillHisto1F(selLabels_[s] + "_RecoPFMET__PseudoCaloMETCut", RecoPFMET[0], weightOnCalib);
                if (RecoPFMET[0] > 170) vcp[s].FillHisto1F(selLabels_[s] + "_PseudoCaloMET__RecoPFMETCut", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_L1MET", L1MET[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_HLTCaloMET", HLTCaloMET[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_HLTCaloMHT", HLTCaloMHT[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_HLTPFMHT", HLTPFMHT[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_HLTPFMET", HLTPFMET[0], weightOnCalib);

                if (*HLT_PFMET120_PFMHT120_IDTight) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    if (PseudoCaloMET[0] > 170) vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET__PseudoCaloMETCut", RecoPFMET[0], weightOnCalib);
                    if (RecoPFMET[0] > 170) vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET__RecoPFMETCut", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if (*HLT_PFHT500_PFMET100_PFMHT100_IDTight) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    if (PseudoCaloMET[0] > 170) vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET__PseudoCaloMETCut", RecoPFMET[0], weightOnCalib);
                    if (RecoPFMET[0] > 170) vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET__RecoPFMETCut", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    if (PseudoCaloMET[0] > 170) vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET__PseudoCaloMETCut", RecoPFMET[0], weightOnCalib);
                    if (RecoPFMET[0] > 170) vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET__RecoPFMETCut", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if (*HLT_MET105_IsoTrk50) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    if (PseudoCaloMET[0] > 170) vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___RecoPFMET__PseudoCaloMETCut", RecoPFMET[0], weightOnCalib);
                    if (RecoPFMET[0] > 170) vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PseudoCaloMET__RecoPFMETCut", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if (*HLT_PFMET120_PFMHT120_IDTight || *HLT_PFHT500_PFMET100_PFMHT100_IDTight
                    || *HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 || *HLT_MET105_IsoTrk50) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMETtrg___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMETtrg___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    if (PseudoCaloMET[0] > 170) vcp[s].FillHisto1F(selLabels_[s] + "_if___orMETtrg___RecoPFMET__PseudoCaloMETCut", RecoPFMET[0], weightOnCalib);
                    if (RecoPFMET[0] > 170) vcp[s].FillHisto1F(selLabels_[s] + "_if___orMETtrg___PseudoCaloMET__RecoPFMETCut", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                
            }

        }
    }


    //  ------------------------------------------------------------------------------------------------------------------------------------------------------------
    //  ------------------------------------------------------------------------------------------------------------------------------------------------------------
    //  ------------------------------------------------------------------------------------------------------------------------------------------------------------
    //  ------------------------------------------------------------------------------------------------------------------------------------------------------------

    for(unsigned int s=0;s<selections_.size();s++) {
            // Trigger Eff vs PUppiMET and vs PseudoMET
        if (selLabels_[s] == "TriggerEffCalib" || selLabels_[s] == "TriggerEffCalib_PseudoMETrescaled") {
            bool BasicSel = false;
            unsigned int hasPassedMuon = 0;
            unsigned int hadPassedHSCP = 0;

            if (*HLT_IsoMu27 && muon_pt.GetSize()==1 && Flag_allMETFilters[0]) BasicSel = true;

            for (unsigned int im = 0; im < muon_pt.GetSize(); im++) {
                if (muon_isTight[im] && muon_pt[im]>30 && muon_pfMiniRelIsoAll[im]<0.15) hasPassedMuon++;
            }

            unsigned int i_track = 0;
            std::vector <unsigned int> i_candPartialHSCP;
            for(unsigned int j=0; j<HSCP_hasTrack.GetSize(); j++){    

                if (!HSCP_hasTrack[j]) continue;

                if ( (Pt[i_track] > 50.0) && (Pt_pseudo[i_track] > 50.0) && (std::abs(Eta[i_track]) < 2.4) && (NbPixelHit_noL1[i_track] >= 2) && (FracOfValidHit[i_track] > 0.8) && 
                (NOM_noL1[i_track] >= 10) && (isHighPurityTrack[i_track] == true) && (normChi2[i_track] < 5.0) && (std::abs(dz[i_track]) < 0.1) && (std::abs(dxy[i_track]) < 0.02) && 
                (miniRelIsoAll[i_track] < 0.02) && (EoP[i_track] < 0.3) && (IsoSumPt_dr03[i_track] < 15) && (ptOverptErrptErr[i_track] < 0.0008) && 
                (Fpix[i_track] < 0.9)) { hadPassedHSCP++; i_candPartialHSCP.push_back(i_track); };

                i_track++;
            }
            
            if (BasicSel && hasPassedMuon==1 && hadPassedHSCP>0) {

                for (unsigned int im = 0; im < muon_pt.GetSize(); im++) {
                    if (muon_isTight[im] && muon_pt[im]>30 && muon_pfMiniRelIsoAll[im]<0.15) {
                        vcp[s].FillHisto1F(selLabels_[s] + "_muon_pt", muon_pt[im], *weightPU);
                        vcp[s].FillHisto1F(selLabels_[s] + "_muon_eta", muon_eta[im], *weightPU);
                        vcp[s].FillHisto1F(selLabels_[s] + "_muon_phi", muon_phi[im], *weightPU);
                    }
                }

                // After a fit performed on Nm1 CaloMET data and MC (max_data/max_MC)
                float isRescaled = -1, weightOnCalib = *weightPU;
                if (selLabels_[s] == "TriggerEffCalib_PseudoMETrescaled") isRescaled = 163.451/146.864; // = 1.113
                else isRescaled = 1;

                for (int iHSCP=0; iHSCP<(int)i_candPartialHSCP.size(); iHSCP++) vcp[s].FillHisto1F(selLabels_[s] + "_PFtrackPT", Pt[i_candPartialHSCP[iHSCP]], weightOnCalib);
                vcp[s].FillHisto2F(selLabels_[s] + "_PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                vcp[s].FillHisto2F(selLabels_[s] + "_PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_PUppiMET", RecoPUppiMET[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_PUppiMET__SF", RecoPUppiMET[0], AppliedWeight1Dpuppi);
                vcp[s].FillHisto1F(selLabels_[s] + "_PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_RecoPFMET", RecoPFMET[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_L1MET", L1MET[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_HLTCaloMET", HLTCaloMET[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_HLTCaloMHT", HLTCaloMHT[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_HLTPFMHT", HLTPFMHT[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_HLTPFMET", HLTPFMET[0], weightOnCalib);


                if (*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if (*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if (*HLT_MET105_IsoTrk50) {
                    for (int iHSCP=0; iHSCP<(int)i_candPartialHSCP.size(); iHSCP++) vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PFtrackPT", Pt[i_candPartialHSCP[iHSCP]], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if ((*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) || (*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100)
                    || (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60) || *HLT_MET105_IsoTrk50) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_PUppiMET__postTrigger", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_PUppiMET__postTrigger__SF", RecoPUppiMET[0], AppliedWeight1Dpuppi);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMETtrg___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMETtrg___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMETtrg___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMETtrg___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMETtrg___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMETtrg___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }

                    // 3 among 4 triggers have passed
                if ((*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) || (*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100)
                    || (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60)) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg1___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg1___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg1___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg1___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg1___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg1___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if ((*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) || (*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100)
                     || *HLT_MET105_IsoTrk50) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg2___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg2___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg2___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg2___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg2___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg2___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if ((*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) || (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60)
                     || *HLT_MET105_IsoTrk50) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg3___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg3___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg3___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg3___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg3___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg3___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if ((*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100) || (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60)
                     || *HLT_MET105_IsoTrk50) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg4___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg4___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg4___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg4___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg4___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg4___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }

            }
        }

            // Trigger Eff vs PUppiMET and vs PseudoMET for Signal MC
        if (selLabels_[s] == "TriggerEffCalib__Signal" || selLabels_[s] == "TriggerEffCalib_PseudoMETrescaled__Signal") {
            unsigned int hadPassedHSCP = 0;

            unsigned int i_track = 0;
            std::vector <unsigned int> i_candPartialHSCP;
            for(unsigned int j=0; j<HSCP_hasTrack.GetSize(); j++){    

                if (!HSCP_hasTrack[j]) continue;

                if ( (Flag_allMETFilters[0] == true) && (Pt[i_track] > 50.0) && (Pt_pseudo[i_track] > 50.0) && (std::abs(Eta[i_track]) < 2.4) && (NbPixelHit_noL1[i_track] >= 2) && 
                (FracOfValidHit[i_track] > 0.8) && (NOM_noL1[i_track] >= 10) && (isHighPurityTrack[i_track] == true) && (normChi2[i_track] < 5.0) && (std::abs(dz[i_track]) < 0.1) &&
                (std::abs(dxy[i_track]) < 0.02) && (miniRelIsoAll[i_track] < 0.02) && (EoP[i_track] < 0.3) && (IsoSumPt_dr03[i_track] < 15) && (ptOverptErrptErr[i_track] < 0.0008) && 
                (Fpix[i_track] > 0.3) && (ptOverptErr[i_track] < 1) && (Ih_Strip[i_track] > 2.9784)) { hadPassedHSCP++; i_candPartialHSCP.push_back(i_track); };

                i_track++;
            }
            
            if (hadPassedHSCP>0) {

                // After a fit performed on Nm1 CaloMET data and MC (max_data/max_MC)
                float isRescaled = -1, weightOnCalib = *weightPU;
                if (selLabels_[s] == "TriggerEffCalib_PseudoMETrescaled__Signal") isRescaled = 163.451/146.864; // = 1.113
                else isRescaled = 1;

                for (int iHSCP=0; iHSCP<(int)i_candPartialHSCP.size(); iHSCP++) vcp[s].FillHisto1F(selLabels_[s] + "_PFtrackPT", Pt[i_candPartialHSCP[iHSCP]], weightOnCalib);
                vcp[s].FillHisto2F(selLabels_[s] + "_PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                vcp[s].FillHisto2F(selLabels_[s] + "_PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_PUppiMET", RecoPUppiMET[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_RecoPFMET", RecoPFMET[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_L1MET", L1MET[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_HLTCaloMET", HLTCaloMET[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_HLTCaloMHT", HLTCaloMHT[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_HLTPFMHT", HLTPFMHT[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_HLTPFMET", HLTPFMET[0], weightOnCalib);

                if (*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if (*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if (*HLT_MET105_IsoTrk50) {
                    for (int iHSCP=0; iHSCP<(int)i_candPartialHSCP.size(); iHSCP++) vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PFtrackPT", Pt[i_candPartialHSCP[iHSCP]], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if ((*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) || (*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100)
                    || (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60) || *HLT_MET105_IsoTrk50) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMETtrg___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMETtrg___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMETtrg___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMETtrg___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMETtrg___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMETtrg___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }

                    // 3 among 4 triggers have passed
                if ((*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) || (*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100)
                    || (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60)) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg1___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg1___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg1___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg1___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg1___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg1___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if ((*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) || (*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100)
                     || *HLT_MET105_IsoTrk50) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg2___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg2___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg2___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg2___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg2___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg2___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if ((*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) || (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60)
                     || *HLT_MET105_IsoTrk50) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg3___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg3___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg3___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg3___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg3___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg3___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if ((*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100) || (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60)
                     || *HLT_MET105_IsoTrk50) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg4___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg4___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg4___PUppiMET", RecoPUppiMET[0], weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMET3a4trg4___PUppiMETNoMu", PUppiMET_NoMu, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg4___PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto2F(selLabels_[s] + "_if___orMET3a4trg4___PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }

            }
        }
    }


    //  ------------------------------------------------------------------------------------------------------------------------------------------------------------
    //  ------------------------------------------------------------------------------------------------------------------------------------------------------------
    //  ------------------------------------------------------------------------------------------------------------------------------------------------------------
    //  ------------------------------------------------------------------------------------------------------------------------------------------------------------




    // using Wjets and Muon dataset
    for(unsigned int s=0;s<selections_.size();s++) {
        if (selLabels_[s] == "CalibPseudoMET_MuWay" || selLabels_[s] == "CalibPseudoMET_MuWay_isRescaled" || selLabels_[s] == "CalibPseudoMET_MuWay_OTHERisRescaled") {
            bool BasicSel = false;
            unsigned int hasPassedMuon = 0;
            unsigned int hadPassedHSCP = 0;

            if (*HLT_IsoMu27 && muon_pt.GetSize()==1 && Flag_allMETFilters[0]) BasicSel = true;

            for (unsigned int im = 0; im < muon_pt.GetSize(); im++) {
                if (muon_isTight[im] && muon_pt[im]>30 && muon_pfMiniRelIsoAll[im]<0.15) hasPassedMuon++;
            }

            unsigned int i_track = 0;
            for(unsigned int j=0; j<HSCP_hasTrack.GetSize(); j++){    

                if (!HSCP_hasTrack[j]) continue;

                if ( (Pt[i_track] > 50.0) && (Pt_pseudo[i_track] > 50.0) && (std::abs(Eta[i_track]) < 2.4) && (NbPixelHit_noL1[i_track] >= 2) && (FracOfValidHit[i_track] > 0.8) && 
                (NOM_noL1[i_track] >= 10) && (isHighPurityTrack[i_track] == true) && (normChi2[i_track] < 5.0) && (std::abs(dz[i_track]) < 0.1) && (std::abs(dxy[i_track]) < 0.02) && 
                (miniRelIsoAll[i_track] < 0.02) && (EoP[i_track] < 0.3) && (IsoSumPt_dr03[i_track] < 15) && (ptOverptErrptErr[i_track] < 0.0008) && 
                (Fpix[i_track] < 0.9)) hadPassedHSCP++;

                i_track++;
            }
            
            if (BasicSel && hasPassedMuon==1 && hadPassedHSCP>0) {

                for (unsigned int im = 0; im < muon_pt.GetSize(); im++) {
                    if (muon_isTight[im] && muon_pt[im]>30 && muon_pfMiniRelIsoAll[im]<0.15) {
                        vcp[s].FillHisto1F(selLabels_[s] + "_muon_pt", muon_pt[im], *weightPU);
                        vcp[s].FillHisto1F(selLabels_[s] + "_muon_eta", muon_eta[im], *weightPU);
                        vcp[s].FillHisto1F(selLabels_[s] + "_muon_phi", muon_phi[im], *weightPU);
                    }
                }

                // After a fit performed on Nm1 CaloMET data and MC (max_data/max_MC)
                float isRescaled = -1, weightOnCalib = *weightPU;
                if (selLabels_[s] == "CalibPseudoMET_MuWay_isRescaled") isRescaled = 163.451/146.864; // = 1.113
                else if (selLabels_[s] == "CalibPseudoMET_MuWay") isRescaled = 1;
                else if (selLabels_[s] == "CalibPseudoMET_MuWay_OTHERisRescaled") isRescaled = (-1.56407e-05*PseudoCaloMET[0]*PseudoCaloMET[0] + 1.01015*PseudoCaloMET[0] + 7.7078)/PseudoCaloMET[0];


                vcp[s].FillHisto1F(selLabels_[s] + "_PUppiMET", RecoPUppiMET[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_RecoPFMET", RecoPFMET[0], weightOnCalib);
                if (PseudoCaloMET[0] > 170) vcp[s].FillHisto1F(selLabels_[s] + "_RecoPFMET__PseudoCaloMETCut", RecoPFMET[0], weightOnCalib);
                if (RecoPFMET[0] > 170) vcp[s].FillHisto1F(selLabels_[s] + "_PseudoCaloMET__RecoPFMETCut", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_L1MET", L1MET[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_HLTCaloMET", HLTCaloMET[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_HLTCaloMHT", HLTCaloMHT[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_HLTPFMHT", HLTPFMHT[0], weightOnCalib);
                vcp[s].FillHisto1F(selLabels_[s] + "_HLTPFMET", HLTPFMET[0], weightOnCalib);

                if (*HLT_PFMET120_PFMHT120_IDTight) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    if (PseudoCaloMET[0] > 170) vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET__PseudoCaloMETCut", RecoPFMET[0], weightOnCalib);
                    if (RecoPFMET[0] > 170) vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET__RecoPFMETCut", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if (*HLT_PFHT500_PFMET100_PFMHT100_IDTight) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    if (PseudoCaloMET[0] > 170) vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET__PseudoCaloMETCut", RecoPFMET[0], weightOnCalib);
                    if (RecoPFMET[0] > 170) vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET__RecoPFMETCut", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    if (PseudoCaloMET[0] > 170) vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET__PseudoCaloMETCut", RecoPFMET[0], weightOnCalib);
                    if (RecoPFMET[0] > 170) vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET__RecoPFMETCut", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if (*HLT_MET105_IsoTrk50) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    if (PseudoCaloMET[0] > 170) vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___RecoPFMET__PseudoCaloMETCut", RecoPFMET[0], weightOnCalib);
                    if (RecoPFMET[0] > 170) vcp[s].FillHisto1F(selLabels_[s] + "_if___HLT_MET105_IsoTrk50___PseudoCaloMET__RecoPFMETCut", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }
                if (*HLT_PFMET120_PFMHT120_IDTight || *HLT_PFHT500_PFMET100_PFMHT100_IDTight
                    || *HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 || *HLT_MET105_IsoTrk50) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMETtrg___PseudoCaloMET", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                    vcp[s].FillHisto1F(selLabels_[s] + "_if___orMETtrg___RecoPFMET", RecoPFMET[0], weightOnCalib);
                    if (PseudoCaloMET[0] > 170) vcp[s].FillHisto1F(selLabels_[s] + "_if___orMETtrg___RecoPFMET__PseudoCaloMETCut", RecoPFMET[0], weightOnCalib);
                    if (RecoPFMET[0] > 170) vcp[s].FillHisto1F(selLabels_[s] + "_if___orMETtrg___PseudoCaloMET__RecoPFMETCut", PseudoCaloMET[0]*isRescaled, weightOnCalib);
                }

            }
        }
    }


    // HSCP loop 
    int i_track = 0;
    for(unsigned int j=0; j<HSCP_hasTrack.GetSize(); j++){    // Every candidate

        if (!HSCP_hasTrack[j]) continue;


        vcp_nosel[0].FillHisto1F("ndEdx_StripOnly", NOM_noL1[i_track] - NbPixelHit_noL1[i_track], *weightPU);
        vcp_nosel[0].FillHisto1F("Nosel_Ptpseudo", Pt_pseudo[i_track], *weightPU);
        vcp_nosel[0].FillHisto1F("Nosel_P", Pt_pseudo[i_track]*cosh(Eta[i_track]), *weightPU);
        vcp_nosel[0].FillHisto1F("Nosel_10000oP", 10000./(Pt_pseudo[i_track]*cosh(Eta[i_track])), *weightPU);
        vcp_nosel[0].FillHisto2F("Nosel_10000oP_vs_Eta", 10000./(Pt_pseudo[i_track]*cosh(Eta[i_track])), Eta[i_track], *weightPU);
        vcp_nosel[0].FillHisto1F("Nosel_eta", Eta[i_track], *weightPU);
        vcp_nosel[0].FillHisto1F("Nosel_phi", Phi[i_track], *weightPU);
        vcp_nosel[0].FillHisto1F("Nosel_PFMiniIso", miniRelIsoAll[i_track], *weightPU);
        vcp_nosel[0].FillHisto1F("Nosel_TrkIso", IsoSumPt_dr03[i_track], *weightPU);
        vcp_nosel[0].FillHisto1F("Nosel_EoverP", EoP[i_track], *weightPU);
        vcp_nosel[0].FillHisto1F("Nosel_PtErr_over_PtPt", ptOverptErrptErr[i_track], *weightPU); 
        if (isMC || isSignal) {
            vcp_nosel[0].FillHisto1F("Nosel_Fpix", Fpix[i_track], *weightPU);
            vcp_nosel[0].FillHisto1F("Nosel_Ih", Ih_Strip[i_track], *weightPU);
            vcp_nosel[0].FillHisto1F("Nosel_Ih_rescaled", Ih_Strip[i_track]*3.2168/3.12692, *weightPU);
            vcp_nosel[0].FillHisto1F("Nosel_Ih_oldCorr", Ih_Strip_oldCorr[i_track], *weightPU);
            vcp_nosel[0].FillHisto1F("Nosel_Ih_noSF", Ih_Strip_noSF[i_track], *weightPU);
            vcp_nosel[0].FillHisto2F("Nosel_10000oP_vs_Ih", 10000./(Pt_pseudo[i_track]*cosh(Eta[i_track])), Ih_Strip[i_track], *weightPU);
        }
        if (!isSignal && !isMC && Fpix[i_track] <= 0.9) { //blind for data in the signal region
            vcp_nosel[0].FillHisto1F("Nosel_Fpix", Fpix[i_track], *weightPU);
            vcp_nosel[0].FillHisto1F("Nosel_Ih", Ih_Strip[i_track], *weightPU);
            vcp_nosel[0].FillHisto1F("Nosel_Ih_rescaled", Ih_Strip[i_track]*3.2168/3.12692, *weightPU);
            vcp_nosel[0].FillHisto1F("Nosel_Ih_oldCorr", Ih_Strip_oldCorr[i_track], *weightPU);
            vcp_nosel[0].FillHisto1F("Nosel_Ih_noSF", Ih_Strip_noSF[i_track], *weightPU);
            vcp_nosel[0].FillHisto2F("Nosel_10000oP_vs_Ih", 10000./(Pt_pseudo[i_track]*cosh(Eta[i_track])), Ih_Strip[i_track], *weightPU);
        }
        vcp_nosel[0].FillHisto1F("Nosel_PtErr_over_Pt", ptOverptErr[i_track], *weightPU);
        vcp_nosel[0].FillHisto1F("Nosel_PthatQCD", PthatQCD[i_track], *weightPU);

        vcp_nosel[0].FillHisto2F("Nosel_Fpix_vs_RunNumber", *Run, Fpix[i_track], *weightPU);
        vcp_nosel[0].FillHisto2F("Nosel_Ih_vs_RunNumber", *Run, Ih_Strip[i_track], *weightPU);
        vcp_nosel[0].FillHisto2F("Nosel_Ih_oldCorr_vs_RunNumber", *Run, Ih_Strip_oldCorr[i_track], *weightPU);
        vcp_nosel[0].FillHisto2F("Nosel_NPV_vs_RunNumber", *Run, *PV_npvsGood, *weightPU);

        vcp_nosel[0].FillHisto2F("Nosel_Fpix_vs_IhnearC", Ih_Strip[i_track], Fpix[i_track], *weightPU);

        if (PseudoCaloMET[0] > 0 && RecoPUppiMET[0] > 0) vcp_nosel[0].FillHisto1F("Noselbis_P", Pt_pseudo[i_track]*cosh(Eta[i_track]), *weightPU);


        if (PF_type[i_track]==13) vcp_nosel[0].FillHisto2F("trackPT_vs_trackPseudoTrackPT__PFmuon", Pt[i_track], Pt_pseudo[i_track], *weightPU);
        if (PF_type[i_track]==211) vcp_nosel[0].FillHisto2F("trackPT_vs_trackPseudoTrackPT__PFpion", Pt[i_track], Pt_pseudo[i_track], *weightPU);
        vcp_nosel[0].FillHisto2F("trackPT_vs_trackPseudoTrackPT", Pt[i_track], Pt_pseudo[i_track], *weightPU);
        vcp_nosel[0].FillHisto2F("trackETA_vs_trackPseudoTrackETA", Eta[i_track], Eta_pseudo[i_track], *weightPU);
        vcp_nosel[0].FillHisto2F("trackPHI_vs_trackPseudoTrackPHI", Phi[i_track], Phi_pseudo[i_track], *weightPU);


        if (trigger && Flag_allMETFilters[0] == true ) vcp_nosel[0].FillHisto1F("PostTrigger_P", Pt_pseudo[i_track]*cosh(Eta[i_track]), *weightPU);
        if (triggerONLY && Flag_allMETFilters[0] == true) vcp_nosel[0].FillHisto1F("PostTriggerbis_P", Pt_pseudo[i_track]*cosh(Eta[i_track]), *weightPU);

        if (trigger && Flag_allMETFilters[0] == true && PseudoCaloMET[0] > 250.) {
            vcp_nosel[0].FillHisto1F("PostTrigger_PseudoCaloMETcut_P", Pt_pseudo[i_track]*cosh(Eta[i_track]), *weightPU);

            vcp_nosel[0].FillHisto1F("PostTriggerWsfup_PseudoCaloMETcut_P", Pt_pseudo[i_track]*cosh(Eta[i_track]), AppliedWeight1Dpseudo_SF_Up);
            vcp_nosel[0].FillHisto1F("PostTriggerWsf_PseudoCaloMETcut_P", Pt_pseudo[i_track]*cosh(Eta[i_track]), AppliedWeight1Dpseudo);
            vcp_nosel[0].FillHisto1F("PostTriggerWsfdown_PseudoCaloMETcut_P", Pt_pseudo[i_track]*cosh(Eta[i_track]), AppliedWeight1Dpseudo_SF_Down);
        }
        if (trigger && Flag_allMETFilters[0] == true && RecoPUppiMET[0] > 150.) {
            vcp_nosel[0].FillHisto1F("PostTrigger_PUppiMETcut_P", Pt_pseudo[i_track]*cosh(Eta[i_track]), *weightPU);

            vcp_nosel[0].FillHisto1F("PostTriggerWsfup_PUppiMETcut_P", Pt_pseudo[i_track]*cosh(Eta[i_track]), AppliedWeight1Dpuppi_SF_Up);
            vcp_nosel[0].FillHisto1F("PostTriggerWsf_PUppiMETcut_P", Pt_pseudo[i_track]*cosh(Eta[i_track]), AppliedWeight1Dpuppi);
            vcp_nosel[0].FillHisto1F("PostTriggerWsfdown_PUppiMETcut_P", Pt_pseudo[i_track]*cosh(Eta[i_track]), AppliedWeight1Dpuppi_SF_Down);
        }
        if (trigger && Flag_allMETFilters[0] == true && (PseudoCaloMET[0] > 250. || RecoPUppiMET[0] > 150.)) {
            vcp_nosel[0].FillHisto1F("PostTrigger_PseudoCaloMETcut_PUppiMETcut_P", Pt_pseudo[i_track]*cosh(Eta[i_track]), *weightPU);

            vcp_nosel[0].FillHisto1F("PostTriggerWsfup_PseudoCaloMETcut_PUppiMETcut_P", Pt_pseudo[i_track]*cosh(Eta[i_track]), AppliedWeight_SF_Up);
            vcp_nosel[0].FillHisto1F("PostTriggerWsf_PseudoCaloMETcut_PUppiMETcut_P", Pt_pseudo[i_track]*cosh(Eta[i_track]), AppliedWeight);
            vcp_nosel[0].FillHisto1F("PostTriggerWsfdown_PseudoCaloMETcut_PUppiMETcut_P", Pt_pseudo[i_track]*cosh(Eta[i_track]), AppliedWeight_SF_Down);
        }



        // ABCD method
        if (isMC) {
            vcp_nosel[0].FillHisto2F("Nosel_pT_vs_Fpixel", Fpix[i_track], Pt_pseudo[i_track], *weightPU);
        }
        else if (!isMC && (Pt_pseudo[i_track] < 70 || Fpix[i_track] <= 0.9)) {
            vcp_nosel[0].FillHisto2F("Nosel_pT_vs_Fpixel", Fpix[i_track], Pt_pseudo[i_track], *weightPU);
        }

        // ONLY FOR NO SEL
        i_track++;
        //-----------------
    }

    for(unsigned int s=0;s<selections_.size();s++) { 
        bool (HSCPSelector::*ptr)(int);
        ptr = selections_[s];

        // COMMENT HERE TO KEEP ALL THE CANDIDATES
        int icand = iCand[s];     // most ionising candidate
        if (icand < 0) continue;
        //----------------------------

        // Setup for up/down variations
        bool isOkforPUppiMET_up = false, isOkforPUppiMET_down = false;
        if (selLabels_[s].find("TestPseudoMETCut_TestPUppiMETCut") != std::string::npos) {
            bool isFine = (this->*ptr)(icand);

            isOkforPUppiMET_down = isFine && (PseudoCaloMET[0] > 250. || PUppiMET_down > 150.);
            isOkforPUppiMET_up   = isFine && (PseudoCaloMET[0] > 250. || PUppiMET_up > 150.);
        }
        else {
            bool isFine = (this->*ptr)(icand);

            isOkforPUppiMET_down = isFine && (PUppiMET_down > 150.);
            isOkforPUppiMET_up   = isFine && (PUppiMET_up > 150.);
        }

        if((this->*ptr)(icand) && (trigger || trigger_varUP || trigger_varDOWN) ) {

            if (trigger) {
                vcp[s].FillHisto1F(selLabels_[s] + "_ndEdx_StripOnly", NOM_noL1[icand] - NbPixelHit_noL1[icand], AppliedWeight1Dpuppi);
                vcp[s].FillHisto1F(selLabels_[s] + "_Ptpseudo", Pt_pseudo[icand], AppliedWeight1Dpuppi);
                vcp[s].FillHisto1F(selLabels_[s] + "_pT", Pt_pseudo[icand], AppliedWeight1Dpuppi);
                vcp[s].FillHisto1F(selLabels_[s] + "_10000oP", 10000./(Pt_pseudo[icand]*cosh(Eta[icand])), AppliedWeight1Dpuppi);
                vcp[s].FillHisto2F(selLabels_[s] + "_10000oP_vs_Eta", 10000./(Pt_pseudo[icand]*cosh(Eta[icand])), Eta[icand], AppliedWeight1Dpuppi);
                vcp[s].FillHisto1F(selLabels_[s] + "_eta", Eta[icand], AppliedWeight1Dpuppi);
                vcp[s].FillHisto1F(selLabels_[s] + "_phi", Phi[icand], AppliedWeight1Dpuppi);
                vcp[s].FillHisto1F(selLabels_[s] + "_PFMiniIso", miniRelIsoAll[icand], AppliedWeight1Dpuppi);
                vcp[s].FillHisto1F(selLabels_[s] + "_TrkIso", IsoSumPt_dr03[icand], AppliedWeight1Dpuppi);
                vcp[s].FillHisto1F(selLabels_[s] + "_EoverP", EoP[icand], AppliedWeight1Dpuppi);
                vcp[s].FillHisto1F(selLabels_[s] + "_PtErr_over_PtPt", ptOverptErrptErr[icand], AppliedWeight1Dpuppi); 
                if (isMC || isSignal) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_Fpix", Fpix[icand], AppliedWeight1Dpuppi);
                    vcp[s].FillHisto1F(selLabels_[s] + "_Ih", Ih_Strip[icand], AppliedWeight1Dpuppi);
                    vcp[s].FillHisto1F(selLabels_[s] + "_Ih_rescaled", Ih_Strip[icand]*3.2168/3.12692, AppliedWeight1Dpuppi);
                    vcp[s].FillHisto1F(selLabels_[s] + "_Ih_oldCorr", Ih_Strip_oldCorr[icand], AppliedWeight1Dpuppi);
                    vcp[s].FillHisto1F(selLabels_[s] + "_Ih_noSF", Ih_Strip_noSF[icand], AppliedWeight1Dpuppi);
                    vcp[s].FillHisto2F(selLabels_[s] + "_10000oP_vs_Ih", 10000./(Pt_pseudo[icand]*cosh(Eta[icand])), Ih_Strip[icand], AppliedWeight1Dpuppi);
                }
                if (!isSignal && !isMC && Fpix[icand] <= 0.9) { //blind for data in the signal region
                    vcp[s].FillHisto1F(selLabels_[s] + "_Fpix", Fpix[icand], AppliedWeight1Dpuppi);
                    vcp[s].FillHisto1F(selLabels_[s] + "_Ih", Ih_Strip[icand], AppliedWeight1Dpuppi);
                    vcp[s].FillHisto1F(selLabels_[s] + "_Ih_rescaled", Ih_Strip[icand]*3.2168/3.12692, AppliedWeight1Dpuppi);
                    vcp[s].FillHisto1F(selLabels_[s] + "_Ih_oldCorr", Ih_Strip_oldCorr[icand], AppliedWeight1Dpuppi);
                    vcp[s].FillHisto1F(selLabels_[s] + "_Ih_noSF", Ih_Strip_noSF[icand], AppliedWeight1Dpuppi);
                    vcp[s].FillHisto2F(selLabels_[s] + "_10000oP_vs_Ih", 10000./(Pt_pseudo[icand]*cosh(Eta[icand])), Ih_Strip[icand], AppliedWeight1Dpuppi);
                }
                vcp[s].FillHisto1F(selLabels_[s] + "_PtErr_over_Pt", ptOverptErr[icand], AppliedWeight1Dpuppi);
                vcp[s].FillHisto1F(selLabels_[s] + "_PthatQCD", PthatQCD[icand], AppliedWeight1Dpuppi);
                vcp[s].FillHisto2F(selLabels_[s] + "_Fpix_vs_RunNumber", *Run, Fpix[icand], AppliedWeight1Dpuppi);
                vcp[s].FillHisto2F(selLabels_[s] + "_Ih_vs_RunNumber", *Run, Ih_Strip[icand], AppliedWeight1Dpuppi);
                vcp[s].FillHisto2F(selLabels_[s] + "_Ih_oldCorr_vs_RunNumber", *Run, Ih_Strip_oldCorr[icand], AppliedWeight1Dpuppi);
                vcp[s].FillHisto2F(selLabels_[s] + "_NPV_vs_RunNumber", *Run, *PV_npvsGood, AppliedWeight1Dpuppi);
                vcp[s].FillHisto2F(selLabels_[s] + "_Fpix_vs_IhnearC", Ih_Strip[icand], Fpix[icand], AppliedWeight1Dpuppi);


                if (selLabels_[s] != "TriggerEffCalib"
                    && selLabels_[s] != "TriggerEffCalib_PseudoMETrescaled"
                    && selLabels_[s] != "TriggerEffCalib__Signal"
                    && selLabels_[s] != "TriggerEffCalib_PseudoMETrescaled__Signal") {

                    vcp[s].FillHisto2F(selLabels_[s] + "_PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0], AppliedWeight1Dpuppi);
                    vcp[s].FillHisto2F(selLabels_[s] + "_PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0], AppliedWeight1Dpuppi);
                    vcp[s].FillHisto1F(selLabels_[s] + "_PUppiMET", RecoPUppiMET[0], AppliedWeight1Dpuppi);
                    vcp[s].FillHisto1F(selLabels_[s] + "_PUppiMETNoMu", PUppiMET_NoMu, AppliedWeight1Dpuppi);
                    vcp[s].FillHisto1F(selLabels_[s] + "_PseudoCaloMET", PseudoCaloMET[0], AppliedWeight1Dpuppi);
                    vcp[s].FillHisto1F(selLabels_[s] + "_RecoPFMET", RecoPFMET[0], AppliedWeight1Dpuppi);

                }

     
                if (selLabels_[s].find("TestPseudoMETCut_TestPUppiMETCut") != std::string::npos) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_P", Pt_pseudo[icand]*cosh(Eta[icand]), AppliedWeight);
                    vcp[s].FillHisto1F(selLabels_[s] + "_Psfup", Pt_pseudo[icand]*cosh(Eta[icand]), AppliedWeight_SF_Up);
                    vcp[s].FillHisto1F(selLabels_[s] + "_Psfdown", Pt_pseudo[icand]*cosh(Eta[icand]), AppliedWeight_SF_Down);
                    vcp[s].FillHisto1F(selLabels_[s] + "_P2", Pt_pseudo[icand]*cosh(Eta[icand]), AppliedWeight);
                }
                else if (selLabels_[s].find("TestPUppiMETCut") != std::string::npos) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_P", Pt_pseudo[icand]*cosh(Eta[icand]), AppliedWeight1Dpuppi);
                    vcp[s].FillHisto1F(selLabels_[s] + "_Psfup", Pt_pseudo[icand]*cosh(Eta[icand]), AppliedWeight1Dpuppi_SF_Up);
                    vcp[s].FillHisto1F(selLabels_[s] + "_Psfdown", Pt_pseudo[icand]*cosh(Eta[icand]), AppliedWeight1Dpuppi_SF_Down);
                    vcp[s].FillHisto1F(selLabels_[s] + "_P2", Pt_pseudo[icand]*cosh(Eta[icand]), AppliedWeight);
                }
                else if (selLabels_[s].find("TestPseudoMETCut") != std::string::npos) {
                    vcp[s].FillHisto1F(selLabels_[s] + "_P", Pt_pseudo[icand]*cosh(Eta[icand]), AppliedWeight1Dpseudo);
                    vcp[s].FillHisto1F(selLabels_[s] + "_Psfup", Pt_pseudo[icand]*cosh(Eta[icand]), AppliedWeight1Dpseudo_SF_Up);
                    vcp[s].FillHisto1F(selLabels_[s] + "_Psfdown", Pt_pseudo[icand]*cosh(Eta[icand]), AppliedWeight1Dpseudo_SF_Down);
                    vcp[s].FillHisto1F(selLabels_[s] + "_P2", Pt_pseudo[icand]*cosh(Eta[icand]), AppliedWeight);
                }
                else {
                    vcp[s].FillHisto1F(selLabels_[s] + "_P", Pt_pseudo[icand]*cosh(Eta[icand]), AppliedWeight);
                    vcp[s].FillHisto1F(selLabels_[s] + "_Psfup", Pt_pseudo[icand]*cosh(Eta[icand]), AppliedWeight_SF_Up);
                    vcp[s].FillHisto1F(selLabels_[s] + "_Psfdown", Pt_pseudo[icand]*cosh(Eta[icand]), AppliedWeight_SF_Down);
                    vcp[s].FillHisto1F(selLabels_[s] + "_P2", Pt_pseudo[icand]*cosh(Eta[icand]), AppliedWeight);
                }


                // ABCD method
                if (isMC) {
                    vcp[s].FillHisto2F(selLabels_[s] + "_pT_vs_Fpixel", Fpix[icand], Pt_pseudo[icand], AppliedWeight1Dpuppi);
                }
                else if (!isMC && (Pt_pseudo[icand] < 70 || Fpix[icand] <= 0.9)) {
                    vcp[s].FillHisto2F(selLabels_[s] + "_pT_vs_Fpixel", Fpix[icand], Pt_pseudo[icand], AppliedWeight1Dpuppi);
                }
            }

            double overP = 10000./(Pt_pseudo[icand]*cosh(Eta[icand]));
            double newWeight = AppliedWeight1Dpuppi;
            double massForRegions = GetMass(Pt_pseudo[icand]*cosh(Eta[icand]), Ih_Strip[icand], K, C);

            if(Pt_pseudo[icand] <= ptcut_ && trigger) {
                if( (Fpix[icand] > fpix3) && (Fpix[icand] <= fpix4) ) vmrp_regionA_3f4[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix3) && (Fpix[icand] <= fpix6) ) vmrp_regionA_3f6[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix3) && (Fpix[icand] <= fpix8) ) vmrp_regionA_3f8[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix3) && (Fpix[icand] <= fpix9) ) vmrp_regionA_3f9[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix4) && (Fpix[icand] <= fpix5) ) vmrp_regionA_4f5[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix5) && (Fpix[icand] <= fpix6) ) vmrp_regionA_5f6[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix6) && (Fpix[icand] <= fpix7) ) vmrp_regionA_6f7[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix6) && (Fpix[icand] <= fpix9) ) vmrp_regionA_6f9[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix7) && (Fpix[icand] <= fpix8) ) vmrp_regionA_7f8[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix8) && (Fpix[icand] <= fpix9) ) vmrp_regionA_8f9[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix9) && (Fpix[icand] <= fpix10) ) vmrp_regionA_9f10[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix99) && (Fpix[icand] <= fpix10) ) vmrp_regionA_99f10[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix999) && (Fpix[icand] <= fpix10) ) vmrp_regionA_999f10[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                
                if( (Fpix[icand] > fpix3) && (Fpix[icand] <= fpix4) ) vmrp_regionB_3f4[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix3) && (Fpix[icand] <= fpix6) ) vmrp_regionB_3f6[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix3) && (Fpix[icand] <= fpix8) ) vmrp_regionB_3f8[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix3) && (Fpix[icand] <= fpix9) ) vmrp_regionB_3f9[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix4) && (Fpix[icand] <= fpix5) ) vmrp_regionB_4f5[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix5) && (Fpix[icand] <= fpix6) ) vmrp_regionB_5f6[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix6) && (Fpix[icand] <= fpix7) ) vmrp_regionB_6f7[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix6) && (Fpix[icand] <= fpix9) ) vmrp_regionB_6f9[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix7) && (Fpix[icand] <= fpix8) ) vmrp_regionB_7f8[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix8) && (Fpix[icand] <= fpix9) ) vmrp_regionB_8f9[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix8) && (Fpix[icand] <= fpix10) ) vmrp_regionB_8f10[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix9) && (Fpix[icand] <= fpix10) ) vmrp_regionB_9f10[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix99) && (Fpix[icand] <= fpix10) ) vmrp_regionB_99f10[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix999) && (Fpix[icand] <= fpix10) ) vmrp_regionB_999f10[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
            }
            else if (Pt_pseudo[icand] > ptcut_ && trigger) {
                if( (Fpix[icand] > fpix3) && (Fpix[icand] <= fpix4) ) vmrp_regionC_3f4[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix3) && (Fpix[icand] <= fpix6) ) vmrp_regionC_3f6[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix3) && (Fpix[icand] <= fpix8) ) vmrp_regionC_3f8[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix3) && (Fpix[icand] <= fpix9) ) vmrp_regionC_3f9[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix4) && (Fpix[icand] <= fpix5) ) vmrp_regionC_4f5[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix5) && (Fpix[icand] <= fpix6) ) vmrp_regionC_5f6[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix6) && (Fpix[icand] <= fpix7) ) vmrp_regionC_6f7[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix6) && (Fpix[icand] <= fpix9) ) vmrp_regionC_6f9[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix7) && (Fpix[icand] <= fpix8) ) vmrp_regionC_7f8[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix8) && (Fpix[icand] <= fpix9) ) vmrp_regionC_8f9[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                
                if( (Fpix[icand] > fpix3) && (Fpix[icand] <= fpix4) ) vmrp_regionD_3f4[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix3) && (Fpix[icand] <= fpix8) ) vmrp_regionD_3f8[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix4) && (Fpix[icand] <= fpix5) ) vmrp_regionD_4f5[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix5) && (Fpix[icand] <= fpix6) ) vmrp_regionD_5f6[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix6) && (Fpix[icand] <= fpix7) ) vmrp_regionD_6f7[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix6) && (Fpix[icand] <= fpix9) ) vmrp_regionD_6f9[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix7) && (Fpix[icand] <= fpix8) ) vmrp_regionD_7f8[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix8) && (Fpix[icand] <= fpix9) ) vmrp_regionD_8f9[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix8) && (Fpix[icand] <= fpix10) ) vmrp_regionD_8f10[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix9) && (Fpix[icand] <= fpix10) ) vmrp_regionD_9f10[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix99) && (Fpix[icand] <= fpix10) ) vmrp_regionD_99f10[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);
                if( (Fpix[icand] > fpix999) && (Fpix[icand] <= fpix10) ) vmrp_regionD_999f10[s].fill(Eta[icand], NOM_noL1[icand], overP, Pt_pseudo[icand], Pterr[icand], Ih_Strip[icand], GStrip[icand], massForRegions, *PV_npvsGood, Fpix[icand], newWeight);


                if (!isSignal && Fpix[icand] > fpix8 && massForRegions >= 1500) {
                    // une seule fois avant la boucle
                    cout << "Mass PFid MET Pt_pseudo Pt Eta NbPixHit NOM FOVH HP chi2 dz dxy PFiso TrkIso EoP ptErr/pt2 ptErr/pt Fpix Ih" << endl;
                    cout << massForRegions << " "
                            << PF_type[icand] << " "
                            << RecoPUppiMET[0] << " " 
                            << Pt_pseudo[icand] << " "
                            << Pt[icand] << " "
                            << Eta[icand] << " "
                            << NbPixelHit_noL1[icand] << " "
                            << NOM_noL1[icand] << " "
                            << FracOfValidHit[icand] << " "
                            << isHighPurityTrack[icand] << " "
                            << normChi2[icand] << " "
                            << dz[icand] << " "
                            << dxy[icand] << " "
                            << miniRelIsoAll[icand] << " "
                            << IsoSumPt_dr03[icand] << " "
                            << EoP[icand] << " "
                            << ptOverptErrptErr[icand] << " "
                            << ptOverptErr[icand] << " "
                            << Fpix[icand] << " "
                            << Ih_Strip[icand] << endl;
                }
            }

            if (isMC) { 
                if (Pt_pseudo[icand] > ptcut_) {
                    if (trigger && Fpix[icand] > fpix9 && Fpix[icand] <= fpix10) {
                        vcp[s].FillHisto1F(selLabels_[s] + "_9fp10_SignalMass_nominal", massForRegions, AppliedWeight1Dpuppi);
                        vcp[s].FillHisto1F(selLabels_[s] + "_9fp10_SignalMass_ATLASdata", findMass(Pt_pseudo[icand]*cosh(Eta[icand]), Ih_Strip[icand], "2024data"), AppliedWeight1Dpuppi);
                        vcp[s].FillHisto1F(selLabels_[s] + "_9fp10_SignalMass_ATLASbkg", findMass(Pt_pseudo[icand]*cosh(Eta[icand]), Ih_Strip[icand], "2024bkg"), AppliedWeight1Dpuppi);
                        vcp[s].FillHisto1F(selLabels_[s] + "_9fp10_SignalMass_ATLASglupion", findMass(Pt_pseudo[icand]*cosh(Eta[icand]), Ih_Strip[icand], "2024glupion"), AppliedWeight1Dpuppi);
                        
                            // PU Up/Down
                        vcp[s].FillHisto1F(selLabels_[s] + "_9fp10_SignalMass_PUUp", massForRegions, AppliedWeight_PU_Up);
                        vcp[s].FillHisto1F(selLabels_[s] + "_9fp10_SignalMass_PUDown", massForRegions, AppliedWeight_PU_Down);
                        
                            // SF Up/Down
                        vcp[s].FillHisto1F(selLabels_[s] + "_9fp10_SignalMass_TriggerSFUp", massForRegions, AppliedWeight1Dpuppi_SF_Up);
                        vcp[s].FillHisto1F(selLabels_[s] + "_9fp10_SignalMass_TriggerSFDown", massForRegions, AppliedWeight1Dpuppi_SF_Down);

                            // Mass method Up/Down
                        vcp[s].FillHisto1F(selLabels_[s] + "_9fp10_SignalMass_KUp", GetMass(Pt_pseudo[icand]*cosh(Eta[icand]),Ih_Strip[icand],K+0.02,C), AppliedWeight1Dpuppi);
                        vcp[s].FillHisto1F(selLabels_[s] + "_9fp10_SignalMass_KDown", GetMass(Pt_pseudo[icand]*cosh(Eta[icand]),Ih_Strip[icand],K-0.02,C), AppliedWeight1Dpuppi);
                        vcp[s].FillHisto1F(selLabels_[s] + "_9fp10_SignalMass_CUp", GetMass(Pt_pseudo[icand]*cosh(Eta[icand]),Ih_Strip[icand],K,C+0.04), AppliedWeight1Dpuppi);
                        vcp[s].FillHisto1F(selLabels_[s] + "_9fp10_SignalMass_CDown", GetMass(Pt_pseudo[icand]*cosh(Eta[icand]),Ih_Strip[icand],K,C-0.04), AppliedWeight1Dpuppi);
                    }
                    
                        // Jet Up/Down
                    if (trigger_varUP && isOkforPUppiMET_up && Fpix[icand] > fpix9 && Fpix[icand] <= fpix10) vcp[s].FillHisto1F(selLabels_[s] + "_9fp10_SignalMass_JetUp", massForRegions, AppliedWeight1Dpuppi);
                    if (trigger_varDOWN && isOkforPUppiMET_down && Fpix[icand] > fpix9 && Fpix[icand] <= fpix10) vcp[s].FillHisto1F(selLabels_[s] + "_9fp10_SignalMass_JetDown", massForRegions, AppliedWeight1Dpuppi);
                }

                if (Pt_pseudo[icand] > ptcut_) {
                    if (trigger && Fpix[icand] > fpix99 && Fpix[icand] <= fpix10) {
                        vcp[s].FillHisto1F(selLabels_[s] + "_99fp10_SignalMass_nominal", massForRegions, AppliedWeight1Dpuppi);
                        vcp[s].FillHisto1F(selLabels_[s] + "_99fp10_SignalMass_ATLASbkg", findMass(Pt_pseudo[icand]*cosh(Eta[icand]), Ih_Strip[icand], "2024data"), AppliedWeight1Dpuppi);
                        vcp[s].FillHisto1F(selLabels_[s] + "_99fp10_SignalMass_ATLASglupion", findMass(Pt_pseudo[icand]*cosh(Eta[icand]), Ih_Strip[icand], "2024glupion"), AppliedWeight1Dpuppi);
                        
                            // PU Up/Down
                        vcp[s].FillHisto1F(selLabels_[s] + "_99fp10_SignalMass_PUUp", massForRegions, AppliedWeight_PU_Up);
                        vcp[s].FillHisto1F(selLabels_[s] + "_99fp10_SignalMass_PUDown", massForRegions, AppliedWeight_PU_Down);
                        
                            // SF Up/Down
                        vcp[s].FillHisto1F(selLabels_[s] + "_99fp10_SignalMass_TriggerSFUp", massForRegions, AppliedWeight1Dpuppi_SF_Up);
                        vcp[s].FillHisto1F(selLabels_[s] + "_99fp10_SignalMass_TriggerSFDown", massForRegions, AppliedWeight1Dpuppi_SF_Down);

                            // Mass method Up/Down
                        vcp[s].FillHisto1F(selLabels_[s] + "_99fp10_SignalMass_KUp", GetMass(Pt_pseudo[icand]*cosh(Eta[icand]),Ih_Strip[icand],K+0.02,C), AppliedWeight1Dpuppi);
                        vcp[s].FillHisto1F(selLabels_[s] + "_99fp10_SignalMass_KDown", GetMass(Pt_pseudo[icand]*cosh(Eta[icand]),Ih_Strip[icand],K-0.02,C), AppliedWeight1Dpuppi);
                        vcp[s].FillHisto1F(selLabels_[s] + "_99fp10_SignalMass_CUp", GetMass(Pt_pseudo[icand]*cosh(Eta[icand]),Ih_Strip[icand],K,C+0.04), AppliedWeight1Dpuppi);
                        vcp[s].FillHisto1F(selLabels_[s] + "_99fp10_SignalMass_CDown", GetMass(Pt_pseudo[icand]*cosh(Eta[icand]),Ih_Strip[icand],K,C-0.04), AppliedWeight1Dpuppi);
                    }

                        // Jet Up/Down
                    if (trigger_varUP && isOkforPUppiMET_up && Fpix[icand] > fpix99 && Fpix[icand] <= fpix10) vcp[s].FillHisto1F(selLabels_[s] + "_99fp10_SignalMass_JetUp", massForRegions, AppliedWeight1Dpuppi);
                    if (trigger_varDOWN && isOkforPUppiMET_down && Fpix[icand] > fpix99 && Fpix[icand] <= fpix10) vcp[s].FillHisto1F(selLabels_[s] + "_99fp10_SignalMass_JetDown", massForRegions, AppliedWeight1Dpuppi);
                }
                if (Pt_pseudo[icand] > ptcut_) {
                    if (trigger && Fpix[icand] > fpix8 && Fpix[icand] <= fpix9) {
                        vcp[s].FillHisto1F(selLabels_[s] + "_8fp9_SignalMass_nominal", massForRegions, AppliedWeight1Dpuppi);
                        vcp[s].FillHisto1F(selLabels_[s] + "_8fp9_SignalMass_ATLASbkg", findMass(Pt_pseudo[icand]*cosh(Eta[icand]), Ih_Strip[icand], "2024data"), AppliedWeight1Dpuppi);
                        vcp[s].FillHisto1F(selLabels_[s] + "_8fp9_SignalMass_ATLASglupion", findMass(Pt_pseudo[icand]*cosh(Eta[icand]), Ih_Strip[icand], "2024glupion"), AppliedWeight1Dpuppi);
                        
                            // PU Up/Down
                        vcp[s].FillHisto1F(selLabels_[s] + "_8fp9_SignalMass_PUUp", massForRegions, AppliedWeight_PU_Up);
                        vcp[s].FillHisto1F(selLabels_[s] + "_8fp9_SignalMass_PUDown", massForRegions, AppliedWeight_PU_Down);
                        
                            // SF Up/Down
                        vcp[s].FillHisto1F(selLabels_[s] + "_8fp9_SignalMass_TriggerSFUp", massForRegions, AppliedWeight1Dpuppi_SF_Up);
                        vcp[s].FillHisto1F(selLabels_[s] + "_8fp9_SignalMass_TriggerSFDown", massForRegions, AppliedWeight1Dpuppi_SF_Down);

                            // Mass method Up/Down
                        vcp[s].FillHisto1F(selLabels_[s] + "_8fp9_SignalMass_KUp", GetMass(Pt_pseudo[icand]*cosh(Eta[icand]),Ih_Strip[icand],K+0.02,C), AppliedWeight1Dpuppi);
                        vcp[s].FillHisto1F(selLabels_[s] + "_8fp9_SignalMass_KDown", GetMass(Pt_pseudo[icand]*cosh(Eta[icand]),Ih_Strip[icand],K-0.02,C), AppliedWeight1Dpuppi);
                        vcp[s].FillHisto1F(selLabels_[s] + "_8fp9_SignalMass_CUp", GetMass(Pt_pseudo[icand]*cosh(Eta[icand]),Ih_Strip[icand],K,C+0.04), AppliedWeight1Dpuppi);
                        vcp[s].FillHisto1F(selLabels_[s] + "_8fp9_SignalMass_CDown", GetMass(Pt_pseudo[icand]*cosh(Eta[icand]),Ih_Strip[icand],K,C-0.04), AppliedWeight1Dpuppi);
                    }

                        // Jet Up/Down
                    if (trigger_varUP && isOkforPUppiMET_up && Fpix[icand] > fpix8 && Fpix[icand] <= fpix9) vcp[s].FillHisto1F(selLabels_[s] + "_8fp9_SignalMass_JetUp", massForRegions, AppliedWeight1Dpuppi);
                    if (trigger_varDOWN && isOkforPUppiMET_down && Fpix[icand] > fpix8 && Fpix[icand] <= fpix9) vcp[s].FillHisto1F(selLabels_[s] + "_8fp9_SignalMass_JetDown", massForRegions, AppliedWeight1Dpuppi);
                }


            }
        
        }
    }


    // Gen part:
    for (unsigned int k = 0; k < GenPart_pt.GetSize(); k++) {
        if (fabs(GenPart_pdgId[k]) > 100000 && fabs(GenPart_charge[k]) > 0) {
            vcp_nosel[0].FillHisto1F("Nosel_GenHSCPcharged_Eta", GenPart_eta[k], *weightPU);
        }
    }


    // --- Pré-calcul de la nature de la paire de gluinos produite ---
    // pairType: 1 = chargé-chargé, 2 = neutre-chargé, 3 = neutre-neutre
    int nChargedHSCP = 0;
    int nHSCP = 0;
    for (unsigned int k = 0; k < GenPart_pt.GetSize(); k++) {
        if (fabs(GenPart_pdgId[k]) > 100000) {
            nHSCP++;
            if (fabs(GenPart_charge[k]) > 0) nChargedHSCP++;
        }
    }
    int pairType = 0;
    if (nHSCP == 2) {
        if      (nChargedHSCP == 2) pairType = 1; // chargé-chargé
        else if (nChargedHSCP == 1) pairType = 2; // neutre-chargé
        else                        pairType = 3; // neutre-neutre
    }


    // trigger eff vs beta + beta dists:
    if (nHSCP==2) {
        vcp_nosel[0].FillHisto1F("Nosel_GenHSCP_beta", (GenPart_beta[0] + GenPart_beta[1])/2, *weightPU);
        if (*HLT_Mu50) vcp_nosel[0].FillHisto1F("Nosel_GenHSCP_beta_ifHLTMu50", (GenPart_beta[0] + GenPart_beta[1])/2, *weightPU);
    }




    // Raw
    vcp_nosel[0].FillHisto1F("Nosel_Gen__PairType", pairType, *weightPU);
    vcp_nosel[0].FillHisto1F("Nosel_Gen__nHSCP", nHSCP, *weightPU);

    // Post trigger
    if (*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) {
        vcp_nosel[0].FillHisto1F("Nosel_Gen__PairType__if_HLT_PFMET120_PFMHT120_IDTight", pairType, *weightPU);
    }
    if (*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100) {
        vcp_nosel[0].FillHisto1F("Nosel_Gen__PairType__if_HLT_PFHT500_PFMET100_PFMHT100_IDTight", pairType, *weightPU);
    }
    if (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60) {
        vcp_nosel[0].FillHisto1F("Nosel_Gen__PairType__if_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", pairType, *weightPU);
    }
    if (*HLT_MET105_IsoTrk50) {
        vcp_nosel[0].FillHisto1F("Nosel_Gen__PairType__if_HLT_MET105_IsoTrk50", pairType, *weightPU);
    }
    if (trigger) {
        vcp_nosel[0].FillHisto1F("Nosel_Gen__PairType__if_ORtrigger", pairType, *weightPU);

        vcp_nosel[0].FillHisto1F("PUppiMET__postTrigger", RecoPUppiMET[0], *weightPU);
        vcp_nosel[0].FillHisto1F("PUppiMET__postTrigger__SF", RecoPUppiMET[0], AppliedWeight1Dpuppi);
    }
    if (*HLT_Mu50) {
        vcp_nosel[0].FillHisto1F("Nosel_Gen__PairType__if_HLT_Mu50", pairType, *weightPU);
    }


    i_track = 0;
    int countMatching = 0;
    for(unsigned int j=0; j<HSCP_hasTrack.GetSize(); j++){    // Every candidate

        if (!HSCP_hasTrack[j]) continue;

        for (unsigned int k = 0; k < GenPart_pt.GetSize(); k++) {
            if (fabs(GenPart_pdgId[k]) > 100000) {
                double dRgen = deltaR(Eta[i_track], Phi[i_track], GenPart_eta[k], GenPart_phi[k]);
                if (dRgen < 0.05) {
                    countMatching++;

                    vcp_nosel[0].FillHisto2F("trackPT_vs_trackPseudoTrackPT__HSCPmatched", Pt[i_track], Pt_pseudo[i_track], *weightPU);
                    vcp_nosel[0].FillHisto1F("GenHSCPmatching__PFType", PF_type[i_track], *weightPU);
                    vcp_nosel[0].FillHisto1F("GenHSCPmatching__HSCP_type", HSCP_type[i_track], *weightPU);
                    if (trigger) vcp_nosel[0].FillHisto1F("GenHSCPmatching__PFType__if_ORtrigger", PF_type[i_track], *weightPU);

                }
            }
        }

        i_track++;
    }
    vcp_nosel[0].FillHisto1F("Nosel_Gen__nMatchedTrack", countMatching, *weightPU);

    // Post trigger + 1 matched HSCP
    if (countMatching==1) {
        vcp_nosel[0].FillHisto1F("Nosel_GenHSCPmatching__PairType__1HSCP", pairType, *weightPU);

        if (*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) {
            vcp_nosel[0].FillHisto1F("Nosel_GenHSCPmatching__PairType__1HSCP_if_HLT_PFMET120_PFMHT120_IDTight", pairType, *weightPU);
        }
        if (*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100) {
            vcp_nosel[0].FillHisto1F("Nosel_GenHSCPmatching__PairType__1HSCP_if_HLT_PFHT500_PFMET100_PFMHT100_IDTight", pairType, *weightPU);
        }
        if (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60) {
            vcp_nosel[0].FillHisto1F("Nosel_GenHSCPmatching__PairType__1HSCP_if_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", pairType, *weightPU);
        }
        if (*HLT_MET105_IsoTrk50) {
            vcp_nosel[0].FillHisto1F("Nosel_GenHSCPmatching__PairType__1HSCP_if_HLT_MET105_IsoTrk50", pairType, *weightPU);
        }
        if (trigger) {
            vcp_nosel[0].FillHisto1F("Nosel_GenHSCPmatching__PairType__1HSCP_if_ORtrigger", pairType, *weightPU);
        }
    }
    // Post trigger + 2 matched HSCP
    else if (countMatching==2) {
        vcp_nosel[0].FillHisto1F("Nosel_GenHSCPmatching__PairType__2HSCP", pairType, *weightPU);

        if (*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) {
            vcp_nosel[0].FillHisto1F("Nosel_GenHSCPmatching__PairType__2HSCP_if_HLT_PFMET120_PFMHT120_IDTight", pairType, *weightPU);
        }
        if (*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100) {
            vcp_nosel[0].FillHisto1F("Nosel_GenHSCPmatching__PairType__2HSCP_if_HLT_PFHT500_PFMET100_PFMHT100_IDTight", pairType, *weightPU);
        }
        if (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60) {
            vcp_nosel[0].FillHisto1F("Nosel_GenHSCPmatching__PairType__2HSCP_if_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", pairType, *weightPU);
        }
        if (*HLT_MET105_IsoTrk50) {
            vcp_nosel[0].FillHisto1F("Nosel_GenHSCPmatching__PairType__2HSCP_if_HLT_MET105_IsoTrk50", pairType, *weightPU);
        }
        if (trigger) {
            vcp_nosel[0].FillHisto1F("Nosel_GenHSCPmatching__PairType__2HSCP_if_ORtrigger", pairType, *weightPU);
        }
    }
    else if (countMatching==0) {
        vcp_nosel[0].FillHisto1F("Nosel_GenHSCPmatching__PairType__0HSCP", pairType, *weightPU);

        if (*HLT_PFMET120_PFMHT120_IDTight && PFMHT>120) {
            vcp_nosel[0].FillHisto1F("Nosel_GenHSCPmatching__PairType__0HSCP_if_HLT_PFMET120_PFMHT120_IDTight", pairType, *weightPU);
        }
        if (*HLT_PFHT500_PFMET100_PFMHT100_IDTight && PFHT>500 && PFMHT>100) {
            vcp_nosel[0].FillHisto1F("Nosel_GenHSCPmatching__PairType__0HSCP_if_HLT_PFHT500_PFMET100_PFMHT100_IDTight", pairType, *weightPU);
        }
        if (*HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && PFMHTNoMu>120 && PFHT>60) {
            vcp_nosel[0].FillHisto1F("Nosel_GenHSCPmatching__PairType__0HSCP_if_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", pairType, *weightPU);
        }
        if (*HLT_MET105_IsoTrk50) {
            vcp_nosel[0].FillHisto1F("Nosel_GenHSCPmatching__PairType__0HSCP_if_HLT_MET105_IsoTrk50", pairType, *weightPU);
        }
        if (trigger) {
            vcp_nosel[0].FillHisto1F("Nosel_GenHSCPmatching__PairType__0HSCP_if_ORtrigger", pairType, *weightPU);
        }
    }




    i_track = 0;
    bool hadPassed_PartialHSCP = false; // No trigger + No PseudoMET cut
    for(unsigned int j=0; j<HSCP_hasTrack.GetSize(); j++){    // Every candidate

        if (!HSCP_hasTrack[j]) continue;

        bool hadPassedHSCP = false;
        if ( (*HLT_PFMET120_PFMHT120_IDTight || *HLT_PFHT500_PFMET100_PFMHT100_IDTight || *HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 || *HLT_MET105_IsoTrk50) && (Flag_allMETFilters[0] == true) && (PseudoCaloMET[0] > 170.) && (Pt[i_track] > 50.0) && (Pt_pseudo[i_track] > 50.0) && (std::abs(Eta[i_track]) < 2.4) && (NbPixelHit_noL1[i_track] >= 2) && (FracOfValidHit[i_track] > 0.8) && (NOM_noL1[i_track] >= 10) && (isHighPurityTrack[i_track] == true) && (normChi2[i_track] < 5.0) && (std::abs(dz[i_track]) < 0.1) && (std::abs(dxy[i_track]) < 0.02) && (miniRelIsoAll[i_track] < 0.02) && (EoP[i_track] < 0.3) && (IsoSumPt_dr03[i_track] < 15) && (ptOverptErrptErr[i_track] < 0.0008) && (Fpix[i_track] > 0.3) && (ptOverptErr[i_track] < 1) && (Ih_Strip[i_track] > 2.9784)) hadPassedHSCP = true;

        for (unsigned int k = 0; k < GenPart_pt.GetSize(); k++) {
            if (fabs(GenPart_pdgId[k]) > 100000) {
                double dRgen = deltaR(Eta[i_track], Phi[i_track], GenPart_eta[k], GenPart_phi[k]);
                if (hadPassedHSCP && dRgen < 0.01 && fabs(GenPart_charge[k]) > 0) {
                    vcp_nosel[0].FillHisto1F("HSCPsel_GenHSCPcharged_Eta", GenPart_eta[k], AppliedWeight1Dpuppi);
                    break;
                }
            }
        }

        if ( (Flag_allMETFilters[0] == true) && (Pt[i_track] > 50.0) && (Pt_pseudo[i_track] > 50.0) && (std::abs(Eta[i_track]) < 2.4) && (NbPixelHit_noL1[i_track] >= 2) && (FracOfValidHit[i_track] > 0.8) && (NOM_noL1[i_track] >= 10) && (isHighPurityTrack[i_track] == true) && (normChi2[i_track] < 5.0) && (std::abs(dz[i_track]) < 0.1) && (std::abs(dxy[i_track]) < 0.02) && (miniRelIsoAll[i_track] < 0.02) && (EoP[i_track] < 0.3) && (IsoSumPt_dr03[i_track] < 15) && (ptOverptErrptErr[i_track] < 0.0008) && (Fpix[i_track] > 0.3) && (ptOverptErr[i_track] < 1) && (Ih_Strip[i_track] > 2.9784)) hadPassed_PartialHSCP = true;

        i_track++;
    }




    // Partial selections
    if (hadPassed_PartialHSCP) {
        vcp_nosel[0].FillHisto2F("HSCPPartialsel_PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0], *weightPU);
        vcp_nosel[0].FillHisto2F("HSCPPartialsel_PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0], *weightPU);
    }
            

    if ((*HLT_PFMET120_PFMHT120_IDTight || *HLT_PFHT500_PFMET100_PFMHT100_IDTight || *HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 || *HLT_MET105_IsoTrk50) 
    && (Flag_allMETFilters[0] == true)) {
        vcp_nosel[0].FillHisto1F("PostTrigger_PseudoCaloMET_SF_rescaled", PseudoCaloMET[0]*163.451/146.864, *weightPU);
        vcp_nosel[0].FillHisto1F("PostTrigger_PseudoCaloMET_SF_NOTrescaled", PseudoCaloMET[0], *weightPU);
    }

    if ((*HLT_PFMET120_PFMHT120_IDTight || *HLT_PFHT500_PFMET100_PFMHT100_IDTight || *HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 || *HLT_MET105_IsoTrk50) 
    && (Flag_allMETFilters[0] == true)) {
        vcp_nosel[0].FillHisto1F("PostTrigger_PseudoCaloMET_rescaled", PseudoCaloMET[0]*163.451/146.864, *weightPU);
        vcp_nosel[0].FillHisto1F("PostTrigger_PseudoCaloMET_NOTrescaled", PseudoCaloMET[0], *weightPU);
    }
    
    vcp_nosel[0].FillHisto2F("Nosel_PseudoMET_vs_PFMET", PseudoCaloMET[0], RecoPFMET[0], *weightPU);
    vcp_nosel[0].FillHisto1F("Nosel_PseudoCaloMET", PseudoCaloMET[0], *weightPU);
    vcp_nosel[0].FillHisto2F("Nosel_PUppiMET_VS_PseudoMET", RecoPUppiMET[0], PseudoCaloMET[0], *weightPU);
    vcp_nosel[0].FillHisto2F("Nosel_PUppiMETNoMu_VS_PseudoMET", PUppiMET_NoMu, PseudoCaloMET[0], *weightPU);
                
    if (isMC) vcp_nosel[0].FillHisto1F("Nosel_PseudoCaloMET_rescaled", PseudoCaloMET[0]*163.451/146.864, *weightPU); // rescale = MAX_data/MAX_MC


    vcp_nosel[0].FillHisto2F("Nosel_PseudoMET_JetSize", PseudoCaloMET[0], Jet_pt.GetSize(), *weightPU);
    vcp_nosel[0].FillHisto2F("Nosel_PUppiMET_JetSize", RecoPUppiMET[0], Jet_pt.GetSize(), *weightPU);

    vcp_nosel[0].FillHisto1F("Nosel_PUppiMET", RecoPUppiMET[0], *weightPU);
    vcp_nosel[0].FillHisto1F("Nosel_PUppiMETup", PUppiMET_up, *weightPU);
    vcp_nosel[0].FillHisto1F("Nosel_PUppiMETdown", PUppiMET_down, *weightPU);

    vcp_nosel[0].FillHisto1F("Trigger", *HLT_FilterOR, *weightPU);
    vcp_nosel[0].FillHisto1F("RecoPFMET", RecoPFMET[0], *weightPU);
    vcp_nosel[0].FillHisto1F("nHSCP", *HSCP_n, *weightPU);
    vcp_nosel[0].FillHisto1F("nPVgood", *PV_npvsGood, *weightPU);



    vcp_nosel[0].FillHisto1F("Nosel_PU_up", *weightPU_Up * *weightPU, 1);
    vcp_nosel[0].FillHisto1F("Nosel_PU_down", *weightPU_Down * *weightPU, 1);
    vcp_nosel[0].FillHisto1F("Nosel_PU", *weightPU, 1);
    vcp_nosel[0].FillHisto1F("Nosel_PU_Factor_up", *weightPU_Up, 1);
    vcp_nosel[0].FillHisto1F("Nosel_PU_Factor_down", *weightPU_Down, 1);




    // dEdx on gluino:
    i_track = 0;
    for(unsigned int j=0; j<HSCP_hasTrack.GetSize(); j++){    // Every candidate

        if (!HSCP_hasTrack[j]) continue;

        if (( (NbPixelHit_noL1[i_track] >= 2) && (FracOfValidHit[i_track] > 0.8) && (NOM_noL1[i_track] >= 10) && (isHighPurityTrack[i_track] == true) && (normChi2[i_track] < 5.0) && (std::abs(dz[i_track]) < 0.5) && (std::abs(dxy[i_track]) < 0.5) && (ptOverptErr[i_track] < 0.25) )) {

            if (Ih_Strip[i_track] >= -20./1500*P[i_track] + 20) {
                vcp_nosel[0].FillHisto2F("dEdX0stripVsP", P[i_track], Ih_Strip[i_track]);
                if (fabs(Eta[i_track]) < 1) vcp_nosel[0].FillHisto2F("dEdX0stripVsP__large", P[i_track], Ih_Strip[i_track]);
            }
        }

        i_track++;
    }

    

   return kTRUE;
}


void HSCPSelector::SlaveTerminate() {
    std::cout << std::endl;
    cout << "Object are saved" << endl;
    for(auto obj: vcp) obj.AddToList(fOutput);
    // The no-selection histograms (cutflows, Nm1_*, Nosel_*, ...) do not depend on the selections.
    // When the selections of a sample are split over several jobs, only one job writes them,
    // otherwise they would be counted several times by hadd.
    if (fillNosel_) for(auto obj: vcp_nosel) obj.AddToList(fOutput);
    else std::cout << "No-selection histograms not written (fillNosel = 0)" << std::endl;

    for(auto obj: vmrp_regionA_3f4) obj.addToList(fOutput);
    for(auto obj: vmrp_regionA_3f6) obj.addToList(fOutput);
    for(auto obj: vmrp_regionA_3f8) obj.addToList(fOutput);
    for(auto obj: vmrp_regionA_3f9) obj.addToList(fOutput);
    for(auto obj: vmrp_regionA_4f5) obj.addToList(fOutput);
    for(auto obj: vmrp_regionA_5f6) obj.addToList(fOutput);
    for(auto obj: vmrp_regionA_6f7) obj.addToList(fOutput);
    for(auto obj: vmrp_regionA_6f9) obj.addToList(fOutput);
    for(auto obj: vmrp_regionA_7f8) obj.addToList(fOutput);
    for(auto obj: vmrp_regionA_8f9) obj.addToList(fOutput);
    for(auto obj: vmrp_regionA_9f10) obj.addToList(fOutput);
    for(auto obj: vmrp_regionA_99f10) obj.addToList(fOutput);
    for(auto obj: vmrp_regionA_999f10) obj.addToList(fOutput);

    for(auto obj: vmrp_regionB_3f4) obj.addToList(fOutput);
    for(auto obj: vmrp_regionB_3f6) obj.addToList(fOutput);
    for(auto obj: vmrp_regionB_3f8) obj.addToList(fOutput);
    for(auto obj: vmrp_regionB_3f9) obj.addToList(fOutput);
    for(auto obj: vmrp_regionB_4f5) obj.addToList(fOutput);
    for(auto obj: vmrp_regionB_5f6) obj.addToList(fOutput);
    for(auto obj: vmrp_regionB_6f7) obj.addToList(fOutput);
    for(auto obj: vmrp_regionB_6f9) obj.addToList(fOutput);
    for(auto obj: vmrp_regionB_7f8) obj.addToList(fOutput);
    for(auto obj: vmrp_regionB_8f9) obj.addToList(fOutput);
    for(auto obj: vmrp_regionB_8f10) obj.addToList(fOutput);
    for(auto obj: vmrp_regionB_9f10) obj.addToList(fOutput);
    for(auto obj: vmrp_regionB_99f10) obj.addToList(fOutput);
    for(auto obj: vmrp_regionB_999f10) obj.addToList(fOutput);

    for(auto obj: vmrp_regionC_3f4) obj.addToList(fOutput);
    for(auto obj: vmrp_regionC_3f6) obj.addToList(fOutput);
    for(auto obj: vmrp_regionC_3f8) obj.addToList(fOutput);
    for(auto obj: vmrp_regionC_3f9) obj.addToList(fOutput);
    for(auto obj: vmrp_regionC_4f5) obj.addToList(fOutput);
    for(auto obj: vmrp_regionC_5f6) obj.addToList(fOutput);
    for(auto obj: vmrp_regionC_6f7) obj.addToList(fOutput);
    for(auto obj: vmrp_regionC_6f9) obj.addToList(fOutput);
    for(auto obj: vmrp_regionC_7f8) obj.addToList(fOutput);
    for(auto obj: vmrp_regionC_8f9) obj.addToList(fOutput);

    for(auto obj: vmrp_regionD_3f4) obj.addToList(fOutput);
    for(auto obj: vmrp_regionD_3f8) obj.addToList(fOutput);
    for(auto obj: vmrp_regionD_4f5) obj.addToList(fOutput);
    for(auto obj: vmrp_regionD_5f6) obj.addToList(fOutput);
    for(auto obj: vmrp_regionD_6f7) obj.addToList(fOutput);
    for(auto obj: vmrp_regionD_6f9) obj.addToList(fOutput);
    for(auto obj: vmrp_regionD_7f8) obj.addToList(fOutput);
    for(auto obj: vmrp_regionD_8f9) obj.addToList(fOutput);
    for(auto obj: vmrp_regionD_8f10) obj.addToList(fOutput);
    for(auto obj: vmrp_regionD_9f10) obj.addToList(fOutput);
    for(auto obj: vmrp_regionD_99f10) obj.addToList(fOutput);
    for(auto obj: vmrp_regionD_999f10) obj.addToList(fOutput);

}

void HSCPSelector::Terminate() {
    //create the output file
    std::cout << std::endl;
    std::cout << "Terminate function" << std::endl;
    fout = new TFile(oFile_.c_str(),"RECREATE");

    //loop over all elements stored
    TIter next(fOutput);
    TObject* object = nullptr;
    while ((object = next())){
        const_cast<TObject*>(object)->Write();
    }

    //write and close the output file
    const_cast<TFile*>(fout)->Write();
    fout->Close();
}