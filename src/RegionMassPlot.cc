#include "RegionMassPlot.h"
#include <TProfile.h>
#include <cmath>

RegionMassPlot::RegionMassPlot(std::string suffix,int etabins,int ihbins,int pbins,int massbins,int fpixbins, float C_parameter)
{

    ih_pt=0;
    ias_pt=0;
    ih_ias=0;
    ih_fpix=0;
    eta_fpix=0;
    oP_fpix=0;
    ih_nhits=0;
    ias_nhits=0;
    eta_pt=0;
    eta_1oP=0;
    eta_p=0;
    eta_pterrOpt=0;
    nhits_pt=0;
    eta_nhits=0;
    ih_eta=0;
    ih_p=0;
    ias_p=0;
    pt_pterroverpt=0;
    ias_eta=0;
    mass_eta=0;
    eta_npv=0;
    p_npv=0;
    ih_npv=0;
    mass=0;
    mass_ih=0;
    mass_p=0;
    suffix_ = suffix;

    initHisto(etabins,ihbins,pbins,massbins,fpixbins,C_parameter);

} 

// The histograms are not deleted here on purpose: RegionMassPlot objects are copied by value
// (the copies share the same pointers) and the histograms are handed over to the selector output list.
RegionMassPlot::~RegionMassPlot(){
}

// Function which intializes the histograms with given binnings 
void RegionMassPlot::initHisto(int& etabins,int& ihbins,int& pbins,int& massbins,int& fpixbins, float& C_parameter)
{
    TH1::AddDirectory(false);
    TH2::AddDirectory(false);

    np = pbins;
    plow = 0;
    pup = 200;

    npt = pbins;
    ptlow = 0;
    ptup = 10000; //instead of 4000

    nih = ihbins;
    ihlow = C_parameter;
    ihup = C_parameter + 5;

    nias = ihbins;
    iaslow = 0;
    iasup = 1;

    neta = etabins;
    etalow = -2.4;
    etaup = 2.4;

    nmass = massbins;
    masslow = 0;
    massup = 4000;

    nfpix = fpixbins;
    fpixlow = 0;
    fpixup = 1.;

    std::string suffix = suffix_;
    
    ih_pt = new TH2F(("ih_pt"+suffix).c_str(),";pt [GeV];I_{h} [MeV/cm]",npt,ptlow,ptup,nih,ihlow,ihup); ih_pt->Sumw2(); ih_pt->SetDirectory(nullptr);
    ias_pt = new TH2F(("ias_pt"+suffix).c_str(),";pt [GeV];I_{as}",npt,ptlow,ptup,nias,iaslow,iasup); ias_pt->Sumw2(); ias_pt->SetDirectory(nullptr);
    ih_ias = new TH2F(("ias_ih"+suffix).c_str(),";I_{as};I_{h} [MeV/cm]",nias,iaslow,iasup,nih,ihlow,ihup); ih_ias->Sumw2(); ih_ias->SetDirectory(nullptr);
    ih_fpix = new TH2F(("fpix_ih"+suffix).c_str(),";F^{pix};I_{h} [MeV/cm]",nfpix,fpixlow,fpixup,nih,ihlow,ihup); ih_fpix->Sumw2(); ih_fpix->SetDirectory(nullptr);
    eta_fpix = new TH2F(("fpix_eta"+suffix).c_str(),";F^{pix};#eta",nfpix,fpixlow,fpixup,neta,etalow,etaup); eta_fpix->Sumw2(); eta_fpix->SetDirectory(nullptr);
    oP_fpix = new TH2F(("oP_fpix"+suffix).c_str(),";10^{4}/p [GeV];F^{pix}",np,plow,pup,nfpix,fpixlow,fpixup); oP_fpix->Sumw2(); oP_fpix->SetDirectory(nullptr);
    ih_nhits = new TH2F(("ih_nhits"+suffix).c_str(),";nhits;I_{h} [MeV/cm]",20,0,20,nih,ihlow,ihup); ih_nhits->Sumw2(); ih_nhits->SetDirectory(nullptr);
    ias_nhits = new TH2F(("ias_nhits"+suffix).c_str(),";nhits;I_{as}",20,0,20,nias,iaslow,iasup); ias_nhits->Sumw2(); ias_nhits->SetDirectory(nullptr);
    eta_pt = new TH2F(("eta_pt"+suffix).c_str(),";pt [GeV];#eta",npt,ptlow,ptup,neta,etalow,etaup); eta_pt->Sumw2(); eta_pt->SetDirectory(nullptr);
    eta_1oP = new TH2F(("eta_1oP"+suffix).c_str(),";10^{4}/p [GeV^{-1}];#eta",np,plow,pup,neta,etalow,etaup); eta_1oP->Sumw2(); eta_1oP->SetDirectory(nullptr);
    eta_p = new TH2F(("eta_p"+suffix).c_str(),";p [GeV];#eta",npt,ptlow,ptup,neta,etalow,etaup); eta_p->Sumw2(); eta_p->SetDirectory(nullptr);
    eta_pterrOpt = new TH2F(("eta_pterrOpt"+suffix).c_str(),";#eta;#frac{#sigma_{pT}}{p_{T}}",neta,etalow,etaup,np,0,1); eta_pterrOpt->Sumw2(); eta_pterrOpt->SetDirectory(nullptr);
    nhits_pt = new TH2F(("nhits_pt"+suffix).c_str(),";pt [GeV];nhits",npt,ptlow,ptup,20,0,20); nhits_pt->Sumw2(); nhits_pt->SetDirectory(nullptr);
    eta_nhits = new TH2F(("eta_nhits"+suffix).c_str(),";nhits;#eta",20,0,20,neta,etalow,etaup); eta_nhits->Sumw2(); eta_nhits->SetDirectory(nullptr);
    ih_eta = new TH2F(("ih_eta"+suffix).c_str(),";#eta;I_{h} [MeV/cm]",neta,etalow,etaup,nih,ihlow,ihup); ih_eta->Sumw2(); ih_eta->SetDirectory(nullptr);
    ih_p = new TH2F(("ih_p"+suffix).c_str(),";10^{4}/p [GeV^{-1}];I_{h} [MeV/cm]",np,plow,pup,nih,ihlow,ihup); ih_p->Sumw2(); ih_p->SetDirectory(nullptr);
    ias_p = new TH2F(("ias_p"+suffix).c_str(),";10^{4}/p [GeV^{-1}];I_{as}",np,plow,pup,nias,iaslow,iasup); ias_p->Sumw2(); ias_p->SetDirectory(nullptr);
    pt_pterroverpt = new TH2F(("pt_pterroverpt"+suffix).c_str(),";p_{T} [GeV];#frac{#sigma_{pT}}{p_{T}}",npt,ptlow,ptup,100,0,1); pt_pterroverpt->Sumw2(); pt_pterroverpt->SetDirectory(nullptr);
    ias_eta = new TH2F(("ias_eta"+suffix).c_str(),";#eta;I_{as}",neta,etalow,etaup,nias,iaslow,iasup); ias_eta->Sumw2(); ias_eta->SetDirectory(nullptr);
    mass_eta = new TH2F(("mass_eta"+suffix).c_str(),";#eta;Mass [GeV]",neta,etalow,etaup,nmass,masslow,massup); mass_eta->Sumw2(); mass_eta->SetDirectory(nullptr);
    eta_npv = new TH2F(("eta_npv"+suffix).c_str(),";npv;#eta",100,0,100,neta,etalow,etaup); eta_npv->Sumw2(); eta_npv->SetDirectory(nullptr);
    p_npv = new TH2F(("p_npv"+suffix).c_str(),";npv;p [GeV]",100,0,100,np,plow,pup); p_npv->Sumw2(); p_npv->SetDirectory(nullptr);
    ih_npv = new TH2F(("ih_npv"+suffix).c_str(),";npv;I_{h} [MeV/cm]",100,0,100,nih,ihlow,ihup); ih_npv->Sumw2(); ih_npv->SetDirectory(nullptr);
    mass = new TH1F(("mass"+suffix).c_str(),";Mass [GeV]",nmass,masslow,massup); mass->Sumw2(); mass->SetDirectory(nullptr);
    mass_p = new TH2F(("mass_p"+suffix).c_str(),";Mass [GeV]",np,plow,pup,nmass,masslow,massup); mass_p->Sumw2(); mass_p->SetDirectory(nullptr);
    mass_ih = new TH2F(("mass_ih"+suffix).c_str(),";Mass [GeV]",nih,ihlow,ihup,nmass,masslow,massup); mass_ih->Sumw2(); mass_ih->SetDirectory(nullptr);
   
    mass->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
}

// Function which fills histograms
void RegionMassPlot::fill(float eta, float nhits, float p, float pt, float pterr, float ih, float ias, float m, float npv, float fpix, float w)
{
   ih_pt->Fill(pt,ih,w);
   ias_pt->Fill(pt,ias,w);
   ih_ias->Fill(ias,ih,w);
   ih_fpix->Fill(fpix,ih,w);
   eta_fpix->Fill(fpix,eta,w);
   oP_fpix->Fill(p,fpix,w);
   ih_nhits->Fill(nhits,ih,w);
   ias_nhits->Fill(nhits,ias,w);
   eta_pt->Fill(pt,eta,w);
   eta_1oP->Fill(p,eta,w);
   eta_p->Fill(pt*cosh(eta),eta,w);
   eta_pterrOpt->Fill(eta,pterr/pt,w);
   nhits_pt->Fill(pt,nhits,w);
   eta_nhits->Fill(nhits,eta,w);
   ih_eta->Fill(eta,ih,w);
  
   ih_p->Fill(p,ih,w);
   ias_p->Fill(p,ias,w);
   mass->Fill(m,w);
   mass_p->Fill(p,m,w);
   mass_ih->Fill(ih,m,w);
   mass_eta->Fill(eta,m,w);

   pt_pterroverpt->Fill(pt,pterr/pt,w);
   ias_eta->Fill(eta,ias,w);
   eta_npv->Fill(npv,eta,w);
   p_npv->Fill(npv,p,w);
   ih_npv->Fill(npv,ih,w);
}


std::vector<TH1*> RegionMassPlot::histograms() const
{
    return {ih_pt, ias_pt, ih_ias, ih_fpix, eta_fpix, oP_fpix, ih_nhits, ias_nhits, eta_pt, eta_1oP, eta_p,
            eta_pterrOpt, nhits_pt, eta_nhits, ih_eta, ih_p, ias_p, pt_pterroverpt, ias_eta, eta_npv, p_npv,
            ih_npv, mass, mass_p, mass_ih, mass_eta};
}


void RegionMassPlot::addFrom(const RegionMassPlot& other)
{
    const std::vector<TH1*> mine = histograms(), theirs = other.histograms();
    for (size_t i = 0; i < mine.size(); ++i)
        if (mine[i] && theirs[i]) mine[i]->Add(theirs[i]);
}


void RegionMassPlot::addToList(TList* list)
{
    auto safeAdd = [&](TObject* o)
    {
        if (!o) return;

        const std::string name = o->GetName();
        if (name.find("CalibPseudoMET") != std::string::npos) return;

        if (!list->FindObject(name.c_str())) list->Add(o);
    };

    safeAdd(ih_pt);
    safeAdd(ias_pt);
    safeAdd(ih_ias);
    safeAdd(ih_fpix);
    safeAdd(eta_fpix);
    safeAdd(oP_fpix);
    safeAdd(ih_nhits);
    safeAdd(ias_nhits);
    safeAdd(eta_pt);
    safeAdd(eta_1oP);
    safeAdd(eta_p);
    safeAdd(eta_pterrOpt);
    safeAdd(nhits_pt);
    safeAdd(eta_nhits);
    safeAdd(ih_eta);
    safeAdd(ih_p);
    safeAdd(ias_p);
    safeAdd(pt_pterroverpt);
    safeAdd(ias_eta);
    safeAdd(eta_npv);
    safeAdd(p_npv);
    safeAdd(ih_npv);
    safeAdd(mass);
    safeAdd(mass_p);
    safeAdd(mass_ih);
    safeAdd(mass_eta);

    // ---- projections ----
    auto h1 = ih_pt->ProjectionY(Form("%s_py", ih_pt->GetName()));
    h1->SetDirectory(nullptr); safeAdd(h1);

    auto h2 = ih_p->ProfileX(Form("%s_pfx", ih_p->GetName()));
    h2->SetDirectory(nullptr); safeAdd(h2);

    auto h3 = ias_p->ProfileX(Form("%s_pfx", ias_p->GetName()));
    h3->SetDirectory(nullptr); safeAdd(h3);

    auto h4 = ih_pt->ProfileX(Form("%s_pfx", ih_pt->GetName()));
    h4->SetDirectory(nullptr); safeAdd(h4);

    auto h5 = ias_pt->ProfileX(Form("%s_pfx", ias_pt->GetName()));
    h5->SetDirectory(nullptr); safeAdd(h5);
}

