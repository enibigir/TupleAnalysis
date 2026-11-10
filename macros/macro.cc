// Driver of the analysis: reads the active line of cfg/configFile.txt, builds the TChain
// of the requested dataset and runs HSCPSelector on it.
//
// Run from the macros/ directory:
//    root -l -b -q macro.cc
//
// To add a dataset: add its name and the location of its ntuples in knownInputs() below,
// and a line in cfg/configFile.txt.

#include <TChain.h>
#include <TROOT.h>
#include <TSystem.h>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

const std::string kConfigFile = "../cfg/configFile.txt";               // relative to macros/
const std::string kDefaultProdDir = "/scratch/ui3_1/gcoulon/HSCP_prod/";   // ntuples and their file lists
const char* const kTreeName   = "HSCPMiniAODAnalyzer/Events";

// kDefaultProdDir, unless the environment variable HSCP_PROD_DIR is set
std::string prodDir() {
    const char* env = std::getenv("HSCP_PROD_DIR");
    if (!env || !*env) return kDefaultProdDir;
    std::string dir = env;
    if (dir.back() != '/') dir += '/';
    return dir;
}

// Location of the ntuples of one dataset, relative to prodDir():
// either a single ROOT file, or a text file listing one ROOT file per line.
struct Input {
    std::string path;
    bool        isList;
};
typedef std::map<std::string, Input> InputMap;

// Signal samples: the dataset <namePrefix><mass> is the ROOT file <filePrefix><mass><fileSuffix>
void addSignal(InputMap& inputs, const std::string& namePrefix, const std::vector<int>& masses,
               const std::string& filePrefix, const std::string& fileSuffix) {
    for (int mass : masses) {
        const std::string m = std::to_string(mass);
        inputs[namePrefix + m] = {filePrefix + m + fileSuffix, false};
    }
}

// Datasets split in several parts: <namePrefix><parts[i]> is the list <listPrefix><first + i>.txt
void addLists(InputMap& inputs, const std::string& namePrefix, const std::vector<std::string>& parts,
              const std::string& listPrefix, int first) {
    for (size_t i = 0; i < parts.size(); ++i) {
        inputs[namePrefix + parts[i]] = {listPrefix + std::to_string(first + static_cast<int>(i)) + ".txt", true};
    }
}

// All the datasets that can be requested in cfg/configFile.txt, for a given code version.
InputMap knownInputs(const std::string& version) {
    InputMap inputs;

    // ---------------- Signal ----------------
    // Gluino: the pythia and madgraph productions are each tied to their own code versions
    if (version == "V19p0") {
        addSignal(inputs, "Gluino_Run3_MET_pythia_", {1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600},
                  "SIGNAL/Gluino_Run3_pythia/Par-M-", "_Code" + version + "_merged.root");
    }
    const std::vector<std::string> madgraphVersions = {"V19p6", "V19p7", "V19p8", "V19p9", "V19p10", "V19p11", "V19p12"};
    if (std::find(madgraphVersions.begin(), madgraphVersions.end(), version) != madgraphVersions.end()) {
        addSignal(inputs, "Gluino_Run3_MET_madgraph_", {1100, 1200, 1300, 1400, 1600, 1800, 2000, 2200, 2400, 2600},
                  "SIGNAL/V19p6/HSCP-Gluino_Par-M-", "_merged.root");
    }
    inputs["Gluino_Run2_MET_madgraph_2000"] = {"SIGNAL/Gluino_Run2_madgraph/Gluino_Run2_MET_madgraph_2000.root", false};

    addSignal(inputs, "Stau_Run3_MET_", {247, 308, 432, 557, 651, 745, 871, 1029, 1218, 1409, 1599},
              "SIGNAL/V20p0/HSCP-Pair-Stau_Par-M-", "_merged.root");
    addSignal(inputs, "Stop_Run3_MET_madgraph_", {700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600},
              "SIGNAL/V21p0/HSCP-Stop_Par-M-", "_merged.root");

    // ---------------- Data: one list per era ----------------
    const std::vector<std::string> eras = {"C", "D", "E", "F", "G", "H", "I"};
    addLists(inputs, "JetMET2024", eras, "JetMET2024/V12p31", 0);
    addLists(inputs, "Mu2024",     eras, "Mu2024/V18p10",     0);
    addLists(inputs, "MuonEG2024", eras, "MuonEG2024/V17p4",  0);

    // ---------------- Simulated backgrounds ----------------
    addLists(inputs, "QCD2024_mu_pt",
             {"15to20", "20to30", "30to50", "50to80", "80to120", "120to170", "170to300",
              "300to470", "470to600", "600to800", "800to1000", "1000"},
             "BKG/QCD2024/V16p2", 1);
    addLists(inputs, "Wjets2024_",
             {"1J_pt40to100", "1J_pt100to200", "1J_pt200to400", "1J_pt400to600", "1J_pt600",
              "2J_pt40to100", "2J_pt100to200", "2J_pt200to400", "2J_pt400to600", "2J_pt600"},
             "BKG/Wjets2024/V14p60", 1);
    inputs["WjetMuNu2024"]     = {"BKG/Wjets2024/V14p13.txt",      true};
    inputs["TTbar2024"]        = {"BKG/TTbar2024/V15p9.txt",       true};
    inputs["TTbarSemiLep2024"] = {"BKG/TTbar1L1Nu2024/V22p0.txt",  true};

    // ---------------- Reduced samples for tests ----------------
    inputs["TestMET2024"]   = {"JetMET2024/testJetMET.txt",        true};
    inputs["TestMuon2024"]  = {"Mu2024/testMu2024.txt",            true};
    inputs["TestMuonEG"]    = {"MuonEG2024/testMuonEG.txt",        true};
    inputs["TestTTbar2024"] = {"BKG/TTbar2024/testTTbar.txt",      true};
    inputs["TestWjets"]     = {"BKG/Wjets2024/testWjets.txt",      true};
    inputs["TestWjetsMuNu"] = {"BKG/Wjets2024/testWjetsMuNu.txt",  true};

    return inputs;
}

// Add the ntuples of one dataset to the chain. Returns false if its file list cannot be opened.
bool fillChain(TChain& chain, const Input& input) {
    const std::string path = prodDir() + input.path;
    if (!input.isList) {
        chain.AddFile(path.c_str());
        return true;
    }

    std::ifstream list(path);
    if (!list.is_open()) {
        std::cerr << "Failed to open file: " << path << std::endl;
        return false;
    }
    std::string line;
    while (std::getline(list, line)) {
        if (!line.empty()) chain.AddFile(line.c_str());
    }
    return true;
}

} // namespace


void macro() {
    gSystem->Load("../libTools.so");
    ROOT::EnableImplicitMT(4);

    // ---------------- Configuration ----------------
    // Lines starting with # are ignored. If several lines are active, the last one is used.
    std::ifstream config(kConfigFile);
    if (!config) {
        std::cerr << "Error when opening config file " << kConfigFile << std::endl;
        return;
    }

    double ptcut = 0.;
    int etabins = 0, ihbins = 0, pbins = 0, massbins = 0, fpixbins = 0;
    std::string dataset, version;

    std::cout << std::endl;
    std::cout << "   Reading config file: " << std::endl;
    std::cout << std::endl;
    std::cout << "pT cut - eta bins - ih bins - p bins - mass bins - FPIXbins -  type - version" << std::endl;
    std::string line;
    while (std::getline(config, line)) {
        if (line.empty() || line[0] == '#') continue;

        double pt;
        int eta, ih, p, mass, fpix;
        std::string name, vers;
        std::stringstream ss(line);
        if (!(ss >> pt >> eta >> ih >> p >> mass >> fpix >> name >> vers)) continue;   // blank or incomplete line

        std::cout << line << std::endl;
        ptcut = pt; etabins = eta; ihbins = ih; pbins = p; massbins = mass; fpixbins = fpix;
        dataset = name; version = vers;
    }
    std::cout << std::endl;

    if (dataset.empty()) {
        std::cerr << "No active line in " << kConfigFile << ". Exiting." << std::endl;
        return;
    }

    // ---------------- Input ntuples ----------------
    const InputMap inputs = knownInputs(version);
    const InputMap::const_iterator input = inputs.find(dataset);
    if (input == inputs.end()) {
        std::cerr << "Dataset " << dataset << " not recognized for version " << version << ". Exiting." << std::endl;
        return;
    }

    TChain chain(kTreeName);
    if (!fillChain(chain, input->second)) return;

    // ---------------- Run the selector ----------------
    // Option string decoded in HSCPSelector::Begin and SlaveBegin
    const std::string binning = std::to_string(ptcut) + "," + std::to_string(etabins) + "," + std::to_string(ihbins) + ","
                              + std::to_string(pbins) + "," + std::to_string(massbins) + "," + std::to_string(fpixbins) + ","
                              + dataset + "," + version;

    std::cout << "Running over dataset : " << dataset << std::endl;
    std::cout << "        code version : " << version << std::endl;
    std::cout << "Defining regions A,B,C,D with pT cut = " << ptcut << std::endl;

    std::cout << "Implicit MT enabled: " << ROOT::IsImplicitMTEnabled() << std::endl;

    chain.SetCacheSize(200 * 1024 * 1024); // 200 MB
    chain.AddBranchToCache("*", true);

    chain.Process("HSCPSelector.C+", binning.c_str());
}
