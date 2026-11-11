// Driver of the analysis: builds the TChain of one sample and runs HSCPSelector on it.
//
// Run from the macros/ directory:
//    root -l -b -q 'macro.cc("JetMET2024C","files.txt","out.root",true,70,96,400,1000,400,20)'
//
//    sample     name, passed to the selector (it chooses the triggers and the K/C constants from it)
//    fileList   text file with one ROOT file per line (absolute paths or URLs)
//    outFile    output ROOT file
//    fillNosel  write the histograms made without selection (cutflows, Nosel_*). Set it on one
//               job only when the selections of a sample are split over several jobs.
//    then the pT cut and the eta, ih, p, mass and Fpix binnings.
//
// The samples, their file lists and the binning are in cfg/datasets.json; launch.sh and
// condor/submit.py resolve them and call this macro. Use "./launch.sh <sample> --dry-run"
// to get the command for one sample.
//
// Environment: HSCP_NTHREADS (implicit-MT threads, default 4; 0 or 1 disables it),
//              HSCP_SF_DIR (trigger scale-factor tables, default ../inc).

#include <TChain.h>
#include <TROOT.h>
#include <TSystem.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

namespace {

const char* const kTreeName = "HSCPMiniAODAnalyzer/Events";

// Number of implicit-MT threads: HSCP_NTHREADS, default 4. 0 or 1 disables implicit MT.
int nThreads() {
    const char* env = std::getenv("HSCP_NTHREADS");
    if (!env || !*env) return 4;
    return std::atoi(env);
}

// Add every non-empty, non-comment line of a file list to the chain. Returns the number of files, -1 if unreadable.
int addFileList(TChain& chain, const std::string& listPath) {
    std::ifstream list(listPath);
    if (!list.is_open()) {
        std::cerr << "Failed to open file list: " << listPath << std::endl;
        return -1;
    }
    int n = 0;
    std::string line;
    while (std::getline(list, line)) {
        while (!line.empty() && (line.back() == ' ' || line.back() == '\t' || line.back() == '\r')) line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        chain.AddFile(line.c_str());
        ++n;
    }
    return n;
}

} // namespace


void macro(const char* sample = "", const char* fileList = "", const char* outFile = "", bool fillNosel = true,
           double ptcut = 70., int etabins = 96, int ihbins = 400, int pbins = 1000, int massbins = 400, int fpixbins = 20) {
    if (!*sample || !*fileList || !*outFile) {
        std::cerr << "usage: macro.cc(\"sample\",\"files.txt\",\"out.root\",fillNosel,ptcut,etabins,ihbins,pbins,massbins,fpixbins)\n"
                  << "       (see ./launch.sh <sample> --dry-run)" << std::endl;
        return;
    }

    gSystem->Load("../libTools.so");

    const int threads = nThreads();
    if (threads > 1) ROOT::EnableImplicitMT(threads);

    TChain chain(kTreeName);
    const int nFiles = addFileList(chain, fileList);
    if (nFiles <= 0) {
        std::cerr << "No file in " << fileList << ". Exiting." << std::endl;
        return;
    }

    // Option string, key=value pairs decoded by HSCPSelector::parseOptions (Begin and SlaveBegin)
    const std::string option = "ptcut=" + std::to_string(ptcut) + ",etabins=" + std::to_string(etabins)
                             + ",ihbins=" + std::to_string(ihbins) + ",pbins=" + std::to_string(pbins)
                             + ",massbins=" + std::to_string(massbins) + ",fpixbins=" + std::to_string(fpixbins)
                             + ",sample=" + sample + ",out=" + outFile + ",nosel=" + (fillNosel ? "1" : "0");

    std::cout << "Running over sample  : " << sample << " (" << nFiles << " files from " << fileList << ")" << std::endl;
    std::cout << "        output file  : " << outFile << std::endl;
    std::cout << "        no-selection histograms : " << (fillNosel ? "written" : "skipped") << std::endl;
    std::cout << "Defining regions A,B,C,D with pT cut = " << ptcut << std::endl;
    std::cout << "Implicit MT enabled: " << ROOT::IsImplicitMTEnabled() << " (" << threads << " threads)" << std::endl;

    chain.SetCacheSize(200 * 1024 * 1024); // 200 MB
    chain.AddBranchToCache("*", true);

    chain.Process("HSCPSelector.C+", option.c_str());
}
