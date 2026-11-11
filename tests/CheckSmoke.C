// Assertions on the output of the smoke test (tests/run_smoke.sh).
//   root -l -b -q 'tests/CheckSmoke.C("out.root", nFiles, expectNosel, "METanalysis_TestPUppiMETCut")'
// The last argument is the selection prefix: <prefix>_Eta2p4 must be in the file, <prefix>_Eta1_2p4 too.
// Exits with status 1 on the first failed check.
#include <TFile.h>
#include <TH1.h>
#include <TSystem.h>
#include <TString.h>
#include <iostream>

static int nFail = 0;
static void expectEntries(TFile& f, const char* name, double expected) {
    TH1* h = dynamic_cast<TH1*>(f.Get(name));
    if (!h) { std::cout << "FAIL  missing histogram " << name << std::endl; ++nFail; return; }
    const double got = h->GetEntries();
    std::cout << (got == expected ? "ok    " : "FAIL  ") << name << ": entries " << got << " expected " << expected << std::endl;
    if (got != expected) ++nFail;
}
static void expectAbsent(TFile& f, const char* name) {
    const bool absent = (f.Get(name) == nullptr);
    std::cout << (absent ? "ok    " : "FAIL  ") << name << (absent ? " absent" : " present but should be absent") << std::endl;
    if (!absent) ++nFail;
}

void CheckSmoke(const char* fileName, int nFiles = 1, bool expectNosel = true,
                const char* prefix = "METanalysis_TestPUppiMETCut") {
    TFile f(fileName);
    if (f.IsZombie()) { std::cout << "FAIL  cannot open " << fileName << std::endl; gSystem->Exit(1); }
    const double n = nFiles;   // events 0 and 2 of each file are selected by the Eta2p4 selections
    const TString sel = TString(prefix) + "_Eta2p4";
    expectEntries(f, Form("%s_pT", sel.Data()), 2 * n);
    expectEntries(f, Form("%s_eta", sel.Data()), 2 * n);
    expectEntries(f, Form("mass_regionD_9fp10_%s", sel.Data()), n);     // event 0: pT 120 > 70, Fpix 0.95
    expectEntries(f, Form("mass_regionD_8fp10_%s", sel.Data()), n);
    expectEntries(f, Form("mass_regionC_8fp9_%s", sel.Data()), 0);
    expectEntries(f, Form("mass_regionA_3fp6_%s", sel.Data()), n);      // event 2: pT 60 < 70, Fpix 0.5
    expectEntries(f, Form("mass_regionB_4fp5_%s", sel.Data()), n);
    expectEntries(f, Form("mass_regionA_9fp10_%s", sel.Data()), 0);
    expectEntries(f, Form("%s_Eta1_2p4_pT", prefix), 0);   // eta 0.5 is outside [1, 2.4]
    if (expectNosel) {
        // per file: events 0 and 2 fill bin 0 + 3 event cuts + 16 candidate cuts = 20 entries each,
        // event 1 (no trigger) fills bin 0 only, event 3 (no candidate) fills bins 0 to 3: 45 entries
        expectEntries(f, "EventCutflow", 45 * n);
        expectEntries(f, "Nosel_eta", 6 * n);                     // 2 + 2 + 2 + 0 tracks per file
    } else {
        expectAbsent(f, "EventCutflow");
        expectAbsent(f, "Nosel_eta");
    }
    std::cout << (nFail ? "SMOKE TEST FAILED" : "SMOKE TEST PASSED") << std::endl;
    gSystem->Exit(nFail ? 1 : 0);
}
