#ifndef HSCPSelector_h
#define HSCPSelector_h

#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TEfficiency.h>
#include <TLorentzVector.h>
#include <vector>
#include <iostream>
#include <limits>
#include <fstream>
#include <map>
#include "ROOT/RConfig.hxx"
#include "TObjString.h"

#include "../inc/RegionMassPlot.h"
#include "../inc/CPlots.h"
#include "../inc/ComputeATLASmass.h"

using namespace std;

// Fpix thresholds of the regions A, B, C, D
constexpr float kFpix3 = 0.3f, kFpix4 = 0.4f, kFpix5 = 0.5f, kFpix6 = 0.6f, kFpix7 = 0.7f, kFpix8 = 0.8f,
                kFpix9 = 0.9f, kFpix99 = 0.99f, kFpix999 = 0.999f, kFpix10 = 1.0f;

// One region of the background estimate in one Fpix slice. Regions A and B hold the candidates with
// pT <= ptcut, C and D those with pT > ptcut; the slice <lo>fp<hi> holds fpixLo < Fpix <= fpixHi.
struct RegionSlice {
   const char* suffix;    // start of the histogram suffix, e.g. "_regionA_3fp4"; the selection label follows
   bool belowPtCut;       // true for A and B, false for C and D
   float fpixLo;
   float fpixHi;
};

// The (region, Fpix slice) pairs, in the order of the output file. To add a slice, add a line here.
const RegionSlice kRegionSlices[] = {
   {"_regionA_3fp4",   true,  kFpix3,   kFpix4},
   {"_regionA_3fp6",   true,  kFpix3,   kFpix6},
   {"_regionA_3fp8",   true,  kFpix3,   kFpix8},
   {"_regionA_3fp9",   true,  kFpix3,   kFpix9},
   {"_regionA_4fp5",   true,  kFpix4,   kFpix5},
   {"_regionA_5fp6",   true,  kFpix5,   kFpix6},
   {"_regionA_6fp7",   true,  kFpix6,   kFpix7},
   {"_regionA_6fp9",   true,  kFpix6,   kFpix9},
   {"_regionA_7fp8",   true,  kFpix7,   kFpix8},
   {"_regionA_8fp9",   true,  kFpix8,   kFpix9},
   {"_regionA_9fp10",  true,  kFpix9,   kFpix10},
   {"_regionA_99fp10", true,  kFpix99,  kFpix10},
   {"_regionA_999fp10",true,  kFpix999, kFpix10},

   {"_regionB_3fp4",   true,  kFpix3,   kFpix4},
   {"_regionB_3fp6",   true,  kFpix3,   kFpix6},
   {"_regionB_3fp8",   true,  kFpix3,   kFpix8},
   {"_regionB_3fp9",   true,  kFpix3,   kFpix9},
   {"_regionB_4fp5",   true,  kFpix4,   kFpix5},
   {"_regionB_5fp6",   true,  kFpix5,   kFpix6},
   {"_regionB_6fp7",   true,  kFpix6,   kFpix7},
   {"_regionB_6fp9",   true,  kFpix6,   kFpix9},
   {"_regionB_7fp8",   true,  kFpix7,   kFpix8},
   {"_regionB_8fp9",   true,  kFpix8,   kFpix9},
   {"_regionB_8fp10",  true,  kFpix8,   kFpix10},
   {"_regionB_9fp10",  true,  kFpix9,   kFpix10},
   {"_regionB_99fp10", true,  kFpix99,  kFpix10},
   {"_regionB_999fp10",true,  kFpix999, kFpix10},

   {"_regionC_3fp4",   false, kFpix3,   kFpix4},
   {"_regionC_3fp6",   false, kFpix3,   kFpix6},
   {"_regionC_3fp8",   false, kFpix3,   kFpix8},
   {"_regionC_3fp9",   false, kFpix3,   kFpix9},
   {"_regionC_4fp5",   false, kFpix4,   kFpix5},
   {"_regionC_5fp6",   false, kFpix5,   kFpix6},
   {"_regionC_6fp7",   false, kFpix6,   kFpix7},
   {"_regionC_6fp9",   false, kFpix6,   kFpix9},
   {"_regionC_7fp8",   false, kFpix7,   kFpix8},
   {"_regionC_8fp9",   false, kFpix8,   kFpix9},

   {"_regionD_3fp4",   false, kFpix3,   kFpix4},
   {"_regionD_3fp8",   false, kFpix3,   kFpix8},
   {"_regionD_4fp5",   false, kFpix4,   kFpix5},
   {"_regionD_5fp6",   false, kFpix5,   kFpix6},
   {"_regionD_6fp7",   false, kFpix6,   kFpix7},
   {"_regionD_6fp9",   false, kFpix6,   kFpix9},
   {"_regionD_7fp8",   false, kFpix7,   kFpix8},
   {"_regionD_8fp9",   false, kFpix8,   kFpix9},
   {"_regionD_8fp10",  false, kFpix8,   kFpix10},
   {"_regionD_9fp10",  false, kFpix9,   kFpix10},
   {"_regionD_99fp10", false, kFpix99,  kFpix10},
   {"_regionD_999fp10",false, kFpix999, kFpix10},
};
const size_t kNRegionSlices = sizeof(kRegionSlices) / sizeof(kRegionSlices[0]);

// Cuts on one candidate with a track, evaluated once per event at the top of Process and read by
// every candidate loop. The thresholds are those of the "common" cuts of cfg/selections.json.
struct TrackCuts {
   bool ptCut = false;        // Pt_pseudo > 50 && Pt > 50
   bool etaCut = false;       // |eta| < 2.4
   bool nPix = false;         // NbPixelHit_noL1 >= 2
   bool nom = false;          // NOM_noL1 >= 10
   bool fovh = false;         // FracOfValidHit > 0.8
   bool highPurity = false;
   bool chi2Cut = false;      // normChi2 < 5
   bool dzCut = false;        // |dz| < 0.1
   bool dxyCut = false;       // |dxy| < 0.02
   bool miniIso = false;      // miniRelIsoAll < 0.02
   bool trkIso = false;       // IsoSumPt_dr03 < 15
   bool eop = false;          // EoP < 0.3
   bool ptErr2 = false;       // ptErr/pt^2 < 0.0008
   bool ptErr = false;        // ptErr/pt < 1
   bool fpixLow = false;      // Fpix > 0.3
   bool fpixCalib = false;    // Fpix < 0.9, used by the trigger-efficiency calibration selections
   bool ihAboveC = false;     // Ih > C of the sample
   bool ihAboveDataC = false; // Ih > 2.9784, the data C, as hard-coded in the partial selections (differs from ihAboveC on MC)

   bool quality()   const { return nPix && nom && fovh && highPurity && chi2Cut && dzCut && dxyCut; }
   bool core()      const { return ptCut && etaCut && quality() && miniIso && trkIso && ptErr2; }   // calibration selections
   bool fullDataC() const { return core() && eop && fpixLow && ptErr && ihAboveDataC; }           // partial HSCP selections
};


class HSCPSelector : public TSelector {
public :
   TTreeReader fReader;  //!the tree reader


   //ADD-HSCP-SELECTION
bool PassHSCPpresel_METanalysis_TestPUppiMETCut_Eta2p4(int hscpIndex);

bool PassHSCPpresel_METanalysis_TestPUppiMETCut_Eta1(int hscpIndex);

bool PassHSCPpresel_METanalysis_TestPUppiMETCut_Eta1_2p4(int hscpIndex);

bool PassHSCPpresel_METanalysis_TestPUppiMETCut_EoP_0p1_Eta2p4(int hscpIndex);

bool PassHSCPpresel_METanalysis_TestPUppiMETCut_EoP_0p1_Eta1(int hscpIndex);

bool PassHSCPpresel_METanalysis_TestPUppiMETCut_EoP_0p1_Eta1_2p4(int hscpIndex);

bool PassHSCPpresel_METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5_Eta2p4(int hscpIndex);

bool PassHSCPpresel_METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5_Eta1(int hscpIndex);

bool PassHSCPpresel_METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5_Eta1_2p4(int hscpIndex);

bool PassHSCPpresel_METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5_EoP_0p1_Eta2p4(int hscpIndex);

bool PassHSCPpresel_METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5_EoP_0p1_Eta1(int hscpIndex);

bool PassHSCPpresel_METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5_EoP_0p1_Eta1_2p4(int hscpIndex);


   int etabins_;
   int ihbins_;
   int pbins_;
   int massbins_;
   int fpixbins_;
   double ptcut_;
   std::string sample_;     // name of the sample, chooses the triggers and the K/C constants
   std::string oFile_;
   bool fillNosel_ = true;   // write the histograms made without selection (see macro.cc)
   std::vector<bool (HSCPSelector::*)(int)> selections_;
   std::vector<string> selLabels_;
   std::vector<CPlots> vcp;
   std::vector<CPlots> vcp_nosel;


   // DATA-MC SF for the signal (derived from the trigg. eff.)
   std::vector <float> SF2D_PseudoMETlovalue, SF2D_PseudoMEThivalue, SF2D_PUppiMETlovalue, SF2D_PUppiMEThivalue, SF2D_Down, SF2D, SF2D_Up;
   std::vector <float> SF1D_PseudoMETvalue, SF1Dpseudo_Down, SF1Dpseudo, SF1Dpseudo_Up;
   std::vector <float> SF1D_PUppiMETvalue, SF1Dpuppi_Down, SF1Dpuppi, SF1Dpuppi_Up;


   TFile* fout;


   // Readers to access the data (delete the ones you do not need).

   // TRIGGERS
   TTreeReaderValue<bool> HLT_Mu50 = {fReader, "HLT_Mu50"};
   TTreeReaderValue<bool> HLT_IsoMu27 = {fReader, "HLT_IsoMu27"};
   TTreeReaderValue<bool> HLT_FilterOR = {fReader, "HLT_FilterOR"};

   TTreeReaderValue<bool> HLT_PFMET120_PFMHT120_IDTight = {fReader, "HLT_PFMET120_PFMHT120_IDTight"};
   TTreeReaderValue<bool> HLT_PFHT500_PFMET100_PFMHT100_IDTight = {fReader, "HLT_PFHT500_PFMET100_PFMHT100_IDTight"};
   TTreeReaderValue<bool> HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = {fReader, "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60"};
   TTreeReaderValue<bool> HLT_MET105_IsoTrk50 = {fReader, "HLT_MET105_IsoTrk50"};

   TTreeReaderValue<bool> HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ = {fReader, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"};
   TTreeReaderValue<bool> HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL = {fReader, "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL"};

   TTreeReaderArray<double> L1MET = {fReader, "L1MET"};
   TTreeReaderArray<float> HLTCaloMET = {fReader, "HLTCaloMET"};
   TTreeReaderArray<float> HLTCaloMHT = {fReader, "HLTCaloMHT"};
   TTreeReaderArray<float> HLTPFMHT = {fReader, "HLTPFMHT"};
   TTreeReaderArray<float> HLTPFMET = {fReader, "HLTPFMET"};

   // EVENT INFO
   TTreeReaderValue<uint32_t> Run = {fReader, "run"};
   TTreeReaderValue<uint32_t> PV_npvsGood = {fReader, "PV_npvsGood"};
   TTreeReaderValue<uint32_t> HSCP_n = {fReader, "HSCP_n"};
   TTreeReaderValue<double> weightPU = {fReader, "weightPU"};


   TTreeReaderValue<double> weightPU_Up = {fReader, "weightPU_Up"};
   TTreeReaderValue<double> weightPU_Down = {fReader, "weightPU_Down"};

   // JETS
   TTreeReaderArray<double> Jet_px = {fReader, "puppijet_px"};
   TTreeReaderArray<double> Jet_py = {fReader, "puppijet_py"};
   TTreeReaderArray<double> Jet_pt = {fReader, "puppijet_pt"};
   TTreeReaderArray<bool> Jet_passJetID = {fReader, "puppijet_passJetID"};
   TTreeReaderArray<bool> Jet_IsRejectedbyJVM = {fReader, "puppijet_IsRejectedbyJVM"}; // true=rejected
   TTreeReaderArray<double> Jet_jesUncTotal = {fReader, "puppijet_jesUncTotal"};
   

   // TRACK INFO
   TTreeReaderArray<int> PF_type = {fReader, "IsoTrack_pfType"};
   TTreeReaderArray<int> HSCP_type = {fReader, "HSCP_type"};
   TTreeReaderArray<bool> HSCP_hasTrack = {fReader, "HSCP_hasTrack"};
   TTreeReaderArray<double> P = {fReader, "IsoTrack_p"};
   TTreeReaderArray<double> Pt = {fReader, "IsoTrack_pt"};
   TTreeReaderArray<double> Pterr = {fReader, "IsoTrack_ptError"};
   TTreeReaderArray<double> Eta = {fReader, "IsoTrack_eta"};
   TTreeReaderArray<double> Phi = {fReader, "IsoTrack_phi"};
   TTreeReaderArray<bool> isHighPurityTrack = {fReader, "IsoTrack_isHighPurityTrack"};
   TTreeReaderArray<double> dz = {fReader, "IsoTrack_dz"};
   TTreeReaderArray<double> dxy = {fReader, "IsoTrack_dxy"};
   TTreeReaderArray<double> normChi2 = {fReader, "IsoTrack_normChi2"};
   TTreeReaderArray<double> EoP = {fReader, "IsoTrack_pfEnergyOverP"};
   TTreeReaderArray<double> ptOverptErrptErr = {fReader, "IsoTrack_ptErrOverPt2"};
   TTreeReaderArray<double> ptOverptErr = {fReader, "IsoTrack_ptErrOverPt"};
   
   TTreeReaderArray<float> Ih_Strip = {fReader, "DeDx_IhStrip"};
   TTreeReaderArray<float> Ih_Strip_noSF = {fReader, "DeDx_IhStrip_noSF"};
   TTreeReaderArray<float> Ih_Strip_oldCorr = {fReader, "DeDx_IhStrip_oldCorr"};

   TTreeReaderArray<uint32_t> NbPixelHit_noL1 = {fReader, "DeDx_PixelNoL1NOM"};
   TTreeReaderArray<uint32_t> NOM_noL1 = {fReader, "DeDx_NoL1NOM"};
   TTreeReaderArray<double> FracOfValidHit = {fReader, "IsoTrack_fractionOfValidHits"};
   
   TTreeReaderArray<float> Fpix = {fReader, "DeDx_FiPixelNoL1"};
   TTreeReaderArray<float> GStrip = {fReader, "DeDx_GiStrip"};

   TTreeReaderArray<float> miniRelIsoAll = {fReader, "IsoTrack_pfMiniRelIsoAll"};
   TTreeReaderArray<float> IsoSumPt_dr03 = {fReader, "IsoTrack_IsoSumPt_dr03"};
   
   TTreeReaderArray<bool> Flag_allMETFilters = {fReader, "Flag_allMETFilters"};
   TTreeReaderArray<double> RecoPFMET = {fReader, "RecoPFMET"};
   TTreeReaderArray<double> PseudoCaloMET = {fReader, "PseudoMET_viaCaloJets"};
   TTreeReaderArray<double> RecoPUppiMET = {fReader, "RecoPuppiMET"};
   TTreeReaderArray<double> RecoPUppiMET_phi = {fReader, "RecoPuppiMET_phi"};

   // Pseudo variable (reco in miniAOD)
   TTreeReaderArray<double> Pt_pseudo = {fReader, "IsoTrack_PseudoTrack_pt"};
   TTreeReaderArray<double> Eta_pseudo = {fReader, "IsoTrack_PseudoTrack_eta"};
   TTreeReaderArray<double> Phi_pseudo = {fReader, "IsoTrack_PseudoTrack_phi"};

   
   // GEN INFO
   TTreeReaderArray<double> GenPart_pt = {fReader, "GenPart_pt"};
   TTreeReaderArray<double> GenPart_eta = {fReader, "GenPart_eta"};
   TTreeReaderArray<double> GenPart_phi = {fReader, "GenPart_phi"};
   TTreeReaderArray<double> GenPart_beta = {fReader, "GenPart_beta"};
   TTreeReaderArray<int> GenPart_pdgId = {fReader, "GenPart_pdgId"};
   TTreeReaderArray<int> GenPart_charge = {fReader, "GenPart_charge"};

   TTreeReaderArray<double> PthatQCD = {fReader, "weightGeneratorBinningValues"};


   // muon info
   TTreeReaderArray<double> muon_pt = {fReader, "muon_pt"};
   TTreeReaderArray<double> muon_eta = {fReader, "muon_eta"};
   TTreeReaderArray<double> muon_phi = {fReader, "muon_phi"};
   TTreeReaderArray<bool> muon_isTight = {fReader, "muon_isTight"};
   TTreeReaderArray<float> muon_pfMiniRelIsoAll = {fReader, "muon_pfMiniRelIsoAll"};
   TTreeReaderArray<bool> muon_isPFMuon = {fReader, "muon_isPFMuon"};


   // electron info
   TTreeReaderArray<double> electron_pt = {fReader, "electron_pt"};
   TTreeReaderArray<double> electron_eta = {fReader, "electron_eta"};
   TTreeReaderArray<double> electron_phi = {fReader, "electron_phi"};
   TTreeReaderArray<float> electron_isTight = {fReader, "electron_isTight"};
   TTreeReaderArray<float> electron_pfMiniRelIsoAll = {fReader, "electron_pfMiniRelIsoAll"};


   // BACKGROUND ESTIMATION METHOD: HISTOGRAMS
   // One vector of RegionMassPlot per entry of kRegionSlices, one element per selection:
   // vmrp_[slice][selection]. The histograms are named <quantity><suffix>_<selection label>.
   std::vector<std::vector<RegionMassPlot>> vmrp_;
   
   HSCPSelector()
   {
      fout = 0;
   }

   virtual ~HSCPSelector(){}
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList   *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();
   std::map<std::string, std::string> parseOptions() const;   // key=value pairs of the option string

   ClassDef(HSCPSelector,0);
};

#endif


#ifdef HSCPSelector_cxx


void HSCPSelector::Init(TTree *tree)
{
   fReader.SetTree(tree);
}

Bool_t HSCPSelector::Notify()
{
   return kTRUE;
}

double deltaR(double eta1, double phi1, double eta2, double phi2) {
   double dphi = phi1 - phi2;
   while (dphi >  M_PI) dphi -= 2*M_PI;
   while (dphi <= -M_PI) dphi += 2*M_PI;

   double deta = eta1 - eta2;
   return std::sqrt(deta*deta + dphi*dphi);
}

#endif // #ifdef HSCPSelector_cxx