#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TFile.h"
#include "TH2.h"
#include "TGraph.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <TPaveText.h>
#include <TSystem.h>
#include <cstdarg>
#include <cstring>
#include <map>
#include <string>

// Output directory of the plots, relative to macros/ (the working directory) or absolute.
// python/plots.py sets it from "output_dir" of cfg/plots.json; the default keeps the old place.
TString gOutDir = "../outputDisplay";

// Legend labels set by python/plots.py from cfg/plots.json, by parameter name of the function ("qcd",
// "gluino", ...) and "lumi"; Label(key, text) gives the configured label, else the text written here.
std::map<std::string, std::string> gLabels;

const char *Label(const char *key, const char *fallback) {
    auto it = gLabels.find(key);
    return it != gLabels.end() ? it->second.c_str() : fallback;
}

// Stops the plot with an error message (python/plots.py reports it as exit code 4).
void Abort(const char *msg) {
    std::cerr << "ERROR: " << msg << std::endl;
    gSystem->Exit(4);
}

// An input that may be left out: no file for an empty path (python/plots.py passes "" for null in
// cfg/plots.json); a path that cannot be opened stops the plot.
TFile *OpenOptional(const char *path) {
    if (!path || !*path) return nullptr;
    TFile *f = TFile::Open(path, "READ");
    if (!f || f->IsZombie()) Abort(Form("cannot open %s", path));
    return f;
}

// gOutDir/<name>, name given printf-style; the parent directory is created if needed.
TString OutPath(const char *fmt, ...) {
    char name[1024];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(name, sizeof(name), fmt, ap);
    va_end(ap);
    TString path = gOutDir + "/" + name;
    gSystem->mkdir(gSystem->GetDirName(path), kTRUE);
    return path;
}

TH1* rebinHisto(TH1* h) {
    double xbins[36]={0.,20.,40.,60.,80.,100.,120.,140.,160.,180.,200.,220.,240.,260.,280.,300.,320.,340.,360.,380.,410.,440.,480.,530.,590.,660.,760.,880.,1030.,1210.,1440.,1730.,2000.,2500.,3200.,4000.};
    std::string newname = h->GetName(); 
    newname += "_rebinned";
    TH1* hres = (TH1*) h->Rebin(35,newname.c_str(),xbins);
    return hres;
}

TH2F* TransposeTH2(const TH2F* h_in) {
    int nx = h_in->GetNbinsX();
    int ny = h_in->GetNbinsY();

    const TAxis* ax = h_in->GetXaxis();
    const TAxis* ay = h_in->GetYaxis();

    std::vector<double> new_xbins, new_ybins;

    for (int j = 1; j <= ny+1; ++j) new_xbins.push_back(ay->GetBinLowEdge(j));
    for (int i = 1; i <= nx+1; ++i) new_ybins.push_back(ax->GetBinLowEdge(i));

    TH2F* h_out = new TH2F(
        Form("%s_transposed", h_in->GetName()),
        Form("%s transposed", h_in->GetTitle()),
        ny, new_xbins.data(),
        nx, new_ybins.data()
    );

    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            double content = h_in->GetBinContent(ix, iy);
            double error = h_in->GetBinError(ix, iy);
            h_out->SetBinContent(iy, ix, content);  // swap X <-> Y
            h_out->SetBinError(iy, ix, error);
        }
    }

    return h_out;
}

double GetMinNonZero(const TH1* h) {
    double min = std::numeric_limits<double>::max();

    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        double c = h->GetBinContent(i);
        if (c > 0 && c < min) min = c;
    }

    if (min == std::numeric_limits<double>::max()) return 0.;

    return min;
}

TCanvas* DrawWithRatio(TH1* h1,
                       TH1* h2,
                       TCanvas* c1,
                       std::string CanvasTitle,
                       std::string RatioTitle,  // h1/h2
                       std::string XaxisTitle,
                       float Xmin,
                       float Xmax) {
    TCanvas* c_new = new TCanvas(CanvasTitle.c_str(), CanvasTitle.c_str(), 800, 600);

    // Define pads
    TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0);
    pad1->SetLeftMargin(0.16);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();

    TPad* pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.315);
    pad2->SetLeftMargin(0.16);pad2->SetBottomMargin(0.33);
    pad2->Draw();


    // Draw upper plot: draw h1 then h2 (use OptionDraw for h1, and "E1 same" for h2)
    pad1->cd();
    c1->DrawClonePad();

    // Draw ratio
    TH1F *h1c = (TH1F*)h1->Clone(TString(h1->GetName()) + "_" + CanvasTitle.c_str() + "_h1");
    TH1F *h2c = (TH1F*)h2->Clone(TString(h2->GetName()) + "_" + CanvasTitle.c_str() + "_h2");
    pad2->cd();
    TH1* h_ratio = (TH1*)h1c->Clone(TString("h_ratio_") + CanvasTitle.c_str());
    h_ratio->Divide(h2c);

    h_ratio->SetTitle("");
    h_ratio->GetYaxis()->SetTitle(RatioTitle.c_str());
    h_ratio->GetXaxis()->SetTitle(XaxisTitle.c_str());
    h_ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
    h_ratio->SetMarkerStyle(8);
    h_ratio->GetYaxis()->SetNdivisions(505);
    h_ratio->GetYaxis()->SetTitleSize(0.08);
    h_ratio->GetYaxis()->SetTitleFont(43);
    h_ratio->GetXaxis()->SetTitleFont(43);
    h_ratio->GetYaxis()->SetLabelFont(43);
    h_ratio->GetXaxis()->SetLabelFont(43);
    h_ratio->GetYaxis()->SetTitleSize(24);  // px
    h_ratio->GetXaxis()->SetTitleSize(24);
    h_ratio->GetYaxis()->SetLabelSize(20);
    h_ratio->GetXaxis()->SetLabelSize(20);
    h_ratio->GetYaxis()->SetTitleOffset(1.1);   // en font pixel, ~1.0 est correct
    h_ratio->SetLineColor(kBlack);
    h_ratio->SetMarkerColor(kBlack);
    gPad->SetTickx(0);
    gStyle->SetOptStat(0);
    h_ratio->LabelsOption("v", "X");
    h_ratio->Draw("E0");
    h_ratio->GetXaxis()->SetRangeUser(Xmin, Xmax);

    TLine* line = new TLine(Xmin, 1, Xmax, 1);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->Draw("same");

    c_new->Update();
    cout << "Canvas " << CanvasTitle << " drawn with ratio of h1: " << h1->GetName() << " to h2: " << h2->GetName() << endl;
    return c_new;
}

TCanvas* DrawWithRatio(TH1* h1,
                       TH1* h2,
                       TH1* h3,
                       TCanvas* c1,
                       std::string CanvasTitle,
                       std::string RatioTitle,  // h1/h2
                       std::string XaxisTitle,
                       float Xmin,
                       float Xmax) {
    TCanvas* c_new = new TCanvas(CanvasTitle.c_str(), CanvasTitle.c_str(), 800, 800);

    // Define pads
    TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0);
    pad1->SetBottomMargin(0.03);
    pad1->Draw();

    TPad* pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.3);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.3);
    pad2->Draw();

    // Draw upper plot: draw h1 then h2 (use OptionDraw for h1, and "E1 same" for h2)
    pad1->cd();
    c1->DrawClonePad();

    // Draw ratio
    TH1F *h1c = (TH1F*)h1->Clone(TString(h1->GetName()) + "_" + CanvasTitle.c_str() + "_h1");
    TH1F *h2c = (TH1F*)h2->Clone(TString(h2->GetName()) + "_" + CanvasTitle.c_str() + "_h2");
    TH1F *h3c = (TH1F*)h3->Clone(TString(h3->GetName()) + "_" + CanvasTitle.c_str() + "_h3");
    pad2->cd();
    TH1* h_ratio = (TH1*)h1c->Clone(TString("h_ratio_") + CanvasTitle.c_str());
    TH1* h_ratiobis = (TH1*)h1c->Clone(TString("h_ratio_") + CanvasTitle.c_str());
    h_ratio->Divide(h2c);
    h_ratiobis->Divide(h3c);

    h_ratio->SetTitle("");
    h_ratio->GetYaxis()->SetTitle(RatioTitle.c_str());
    h_ratio->GetXaxis()->SetTitle(XaxisTitle.c_str());
    h_ratio->GetYaxis()->SetRangeUser(0, 2);
    h_ratio->SetMarkerStyle(8);
    h_ratiobis->SetMarkerStyle(22);
    h_ratio->GetYaxis()->SetNdivisions(505);
    h_ratio->GetYaxis()->SetTitleSize(0.08);
    h_ratio->GetYaxis()->SetTitleOffset(0.3);
    h_ratio->GetXaxis()->SetTitleSize(0.08);
    h_ratio->GetXaxis()->SetTitleOffset(1);
    h_ratio->GetYaxis()->SetLabelSize(0.06);
    h_ratio->GetXaxis()->SetLabelSize(0.06);
    h_ratio->SetLineColor(kBlack);
    h_ratio->SetMarkerColor(kBlack);
    h_ratiobis->SetLineColor(kRed);
    h_ratiobis->SetMarkerColor(kRed);
    gPad->SetTickx(0);
    gStyle->SetOptStat(0);
    h_ratio->LabelsOption("v", "X");
    h_ratio->Draw("E0");
    h_ratiobis->Draw("E0 same");
    h_ratio->GetXaxis()->SetRangeUser(Xmin, Xmax);

    // add legend
    TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(h_ratio, "rescaled", "pe");
    leg->AddEntry(h_ratiobis, "no rescaled", "pe");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();

    TLine* line = new TLine(Xmin, 1, Xmax, 1);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->Draw("same");

    c_new->Update();
    cout << "Canvas " << CanvasTitle << " drawn with ratio of h1: " << h1->GetName() << " to h2: " << h2->GetName() << " and h3: " << h3->GetName() << endl;
    return c_new;
}

TCanvas *DrawCanvas(TH1* h,
                    std::string CanvasTitle,
                    std::string XaxisTitle,
                    std::string YaxisTitle,
                    std::string OptionDraw,
                    float sizetitleY,
                    float Xmin,
                    float Xmax,
                    int color,
                    float Ymin = 0, 
                    float Ymax = -1, 
                    bool isLOGy = false) {

    TCanvas* c_new = new TCanvas(CanvasTitle.c_str(), CanvasTitle.c_str(), 800, 600);
    c_new->cd(); 
    c_new->SetLeftMargin(0.16);
    c_new->SetBottomMargin(0.16);

    c_new->SetTitle("");

    gStyle->SetOptStat(0);
    
    TH1* hc = (TH1*)h->Clone(TString(h->GetName()) + "_" + CanvasTitle.c_str());
    hc->GetYaxis()->SetTitleOffset(1.2);
    hc->SetTitle("");
    hc->GetYaxis()->SetTitleSize(sizetitleY);
    hc->GetXaxis()->SetTitleSize(0.06);
    hc->GetXaxis()->SetTitleOffset(0.9);
    hc->GetYaxis()->SetTitleOffset(1);
    hc->GetXaxis()->SetLabelSize(0.05);
    hc->GetYaxis()->SetLabelSize(0.05);

    hc->GetXaxis()->SetTitle(XaxisTitle.c_str());
    hc->GetYaxis()->SetTitle(YaxisTitle.c_str());
    hc->SetLineColor(color);
    hc->SetMarkerColor(color);
    hc->SetMarkerStyle(20);
    hc->Draw(OptionDraw.c_str());
    hc->GetXaxis()->SetRangeUser(Xmin, Xmax);
    if (Ymax == -1) Ymax = 1.2*hc->GetMaximum();
    //if (isLOGy) Ymin = hc->GetMinimum()*0.8;
    hc->GetYaxis()->SetRangeUser(Ymin, Ymax);
    if (isLOGy) c_new->SetLogy();

    c_new->Update();
    cout << "Canvas " << CanvasTitle << " drawn with h: " << h->GetName() << endl;
    return c_new;
}

TCanvas *DrawCanvas(TH2* h,
                    std::string CanvasTitle,
                    std::string XaxisTitle,
                    std::string YaxisTitle,
                    std::string ZaxisTitle,
                    std::string OptionDraw,
                    float Xmin,
                    float Xmax, 
                    float Ymin, 
                    float Ymax,
                    float Zmin = 0,
                    float Zmax = 0,
                    bool isLOGz = false) {

    TCanvas* c_new = new TCanvas(CanvasTitle.c_str(), CanvasTitle.c_str(), 800, 600);
    c_new->cd(); 
    c_new->SetLeftMargin(0.16);
    c_new->SetBottomMargin(0.16);
    c_new->SetRightMargin(0.16);
    

    c_new->SetTitle("");

    gStyle->SetOptStat(0);
    
    TH2* hc = (TH2*)h->Clone(TString(h->GetName()) + "_" + CanvasTitle.c_str());
    hc->GetYaxis()->SetTitleOffset(1.2);
    hc->SetTitle("");
    hc->GetYaxis()->SetTitleSize(0.06);
    hc->GetXaxis()->SetTitleSize(0.06);
    hc->GetZaxis()->SetTitleSize(0.06);
    hc->GetXaxis()->SetTitleOffset(0.9);
    hc->GetYaxis()->SetTitleOffset(1);
    hc->GetXaxis()->SetLabelSize(0.05);
    hc->GetYaxis()->SetLabelSize(0.05);
    hc->GetZaxis()->SetLabelSize(0.05);
    hc->GetZaxis()->SetTitleOffset(0.9);

    hc->GetXaxis()->SetTitle(XaxisTitle.c_str());
    hc->GetYaxis()->SetTitle(YaxisTitle.c_str());
    hc->GetZaxis()->SetTitle(ZaxisTitle.c_str());
    hc->Draw(OptionDraw.c_str());
    gStyle->SetPaintTextFormat("4.2f");
    hc->GetXaxis()->SetRangeUser(Xmin, Xmax);
    hc->GetYaxis()->SetRangeUser(Ymin, Ymax);
    if (Zmin!=0 || Zmax!=0) hc->GetZaxis()->SetRangeUser(Zmin, Zmax);
    if (isLOGz) c_new->SetLogz();

    c_new->Update();
    cout << "Canvas " << CanvasTitle << " drawn with h: " << h->GetName() << endl;
    return c_new;
}

TCanvas *DrawCanvas(TH1* h1,
                    TH1* h2,
                    std::string CanvasTitle,
                    std::string XaxisTitle,
                    std::string YaxisTitle,
                    std::string OptionDraw_h1,
                    std::string OptionDraw_h2,
                    bool isLegend,
                    std::string legtitle_h1,
                    std::string legtitle_h2,
                    std::string legdraw_h1,
                    std::string legdraw_h2,
                    float sizetitleY,
                    float Xmin,
                    float Xmax, 
                    float Ymin = 0, 
                    float Ymax = -1, 
                    bool isLOGy = false,
                    bool isMCfill_h2 = false) {

    TCanvas* c_new = new TCanvas(CanvasTitle.c_str(), CanvasTitle.c_str(), 800, 600);
    c_new->SetLeftMargin(0.16);
    c_new->SetBottomMargin(0.16);
    

    TH1* h1c = (TH1*)h1->Clone(TString(h1->GetName()) + "_" + CanvasTitle.c_str());
    TH1* h2c = (TH1*)h2->Clone(TString(h2->GetName()) + "_" + CanvasTitle.c_str());

    h1c->GetYaxis()->SetTitleOffset(1.0);
    h1c->SetTitle("");
    h1c->GetYaxis()->SetTitleSize(sizetitleY);
    h1c->GetXaxis()->SetTitleSize(0.06);
    h1c->GetXaxis()->SetTitleOffset(0.9);
    h1c->GetYaxis()->SetTitleOffset(1);
    h1c->GetXaxis()->SetLabelSize(0.05);
    h1c->GetYaxis()->SetLabelSize(0.05);
    h1c->GetXaxis()->SetTitle(XaxisTitle.c_str());
    h1c->GetYaxis()->SetTitle(YaxisTitle.c_str());

    h1c->SetLineColor(kRed);
    h1c->SetMarkerColor(kRed);
    h2c->SetLineColor(kBlue);
    h2c->SetMarkerColor(kBlue);

    if (isMCfill_h2) {
        h1c->SetLineColor(kBlack);
        h1c->SetMarkerColor(kBlack);
        h1c->SetMarkerStyle(20);
        h2c->SetLineColor(kBlue-7);
        h2c->SetMarkerColor(kBlue-7);
        h2c->SetFillStyle(1001);
        h2c->SetFillColor(kBlue-7);
    }

    gStyle->SetOptStat(0);

    if (isMCfill_h2) {
        h1c->Draw(OptionDraw_h1.c_str());
        std::string OptionDraw_h1_filled = OptionDraw_h1 + " same";
        h2c->Draw(OptionDraw_h2.c_str());
        h1c->Draw(OptionDraw_h1_filled.c_str());
    }
    else {
        h1c->Draw(OptionDraw_h1.c_str());
        h2c->Draw(OptionDraw_h2.c_str());
    }

    h1c->GetXaxis()->SetRangeUser(Xmin, Xmax);
    if (Ymax == -1) Ymax = 1.2*std::max(h1->GetMaximum(), h2->GetMaximum());
    if (isLOGy) Ymin = 0.8*std::min(GetMinNonZero(h1c), GetMinNonZero(h2c));
    h1c->GetYaxis()->SetRangeUser(Ymin, Ymax);

    TLegend* leg = new TLegend(0.4, 0.2, 0.6, 0.4);
    leg->AddEntry(h1c, legtitle_h1.c_str(), legdraw_h1.c_str());
    leg->AddEntry(h2c, legtitle_h2.c_str(), legdraw_h2.c_str());
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    if (isLegend) leg->Draw();

    if (isLOGy) c_new->SetLogy();

    c_new->Update();
    cout << "Canvas " << CanvasTitle << " drawn with h1: " << h1->GetName() << " and h2: " << h2->GetName() << endl;
    return c_new;
}

TCanvas *DrawCanvas(TH1* h1,
                    TH1* h2,
                    TH1* h3,
                    std::string CanvasTitle,
                    std::string XaxisTitle,
                    std::string YaxisTitle,
                    std::string OptionDraw_h1,
                    std::string OptionDraw_h2,
                    std::string OptionDraw_h3,
                    bool isLegend,
                    std::string legtitle_h1,
                    std::string legtitle_h2,
                    std::string legtitle_h3,
                    std::string legdraw_h1,
                    std::string legdraw_h2,
                    std::string legdraw_h3,
                    float Xmin,
                    float Xmax, 
                    float Ymin = 0, 
                    float Ymax = -1, 
                    bool isLOGy = false,
                    bool isMCfill_h2 = false) {

    TCanvas* c_new = new TCanvas(CanvasTitle.c_str(), CanvasTitle.c_str(), 800, 600);

    TH1* h1c = (TH1*)h1->Clone(TString(h1->GetName()) + "_" + CanvasTitle.c_str());
    TH1* h2c = (TH1*)h2->Clone(TString(h2->GetName()) + "_" + CanvasTitle.c_str());
    TH1* h3c = (TH1*)h3->Clone(TString(h3->GetName()) + "_" + CanvasTitle.c_str());

    h1c->SetTitle(CanvasTitle.c_str());
    h1c->GetXaxis()->SetTitle(XaxisTitle.c_str());
    h1c->GetYaxis()->SetTitle(YaxisTitle.c_str());

    h1c->SetLineColor(kRed);
    h1c->SetMarkerColor(kRed);
    h2c->SetLineColor(kBlue);
    h2c->SetMarkerColor(kBlue);
    h3c->SetLineColor(kGreen+2);
    h3c->SetMarkerColor(kGreen+2);
        
    if (isMCfill_h2) {
        h1c->SetLineColor(kBlack);
        h1c->SetMarkerColor(kBlack);
        h1c->SetMarkerStyle(20);
        h2c->SetLineColor(kBlue-7);
        h2c->SetMarkerColor(kBlue-7);
        h2c->SetFillStyle(1001);
        h2c->SetFillColor(kBlue-7);
        h3c->SetLineColor(kRed);
        h3c->SetMarkerColor(kRed);
    }

    if (isMCfill_h2) {
        h1c->Draw(OptionDraw_h1.c_str());
        std::string OptionDraw_h1_filled = OptionDraw_h1 + " same";
        h2c->Draw(OptionDraw_h2.c_str());
        h1c->Draw(OptionDraw_h1_filled.c_str());
        h3c->Draw(OptionDraw_h3.c_str());
    }
    else {
        h1c->Draw(OptionDraw_h1.c_str());
        h2c->Draw(OptionDraw_h2.c_str());
        h3c->Draw(OptionDraw_h3.c_str());
    }

    h1c->GetXaxis()->SetRangeUser(Xmin, Xmax);
    if (Ymax == -1) Ymax = 1.2*std::max(h1->GetMaximum(), h2->GetMaximum());
    if (isLOGy) Ymin = 0.8*std::min(GetMinNonZero(h1c), GetMinNonZero(h2c));
    h1c->GetYaxis()->SetRangeUser(Ymin, Ymax);

    TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(h1c, legtitle_h1.c_str(), legdraw_h1.c_str());
    leg->AddEntry(h2c, legtitle_h2.c_str(), legdraw_h2.c_str());
    leg->AddEntry(h3c, legtitle_h3.c_str(), legdraw_h3.c_str());
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    if (isLegend) leg->Draw();

    if (isLOGy) c_new->SetLogy();

    c_new->Update();
    cout << "Canvas " << CanvasTitle << " drawn with h1: " << h1->GetName() << " and h2: " << h2->GetName() << " and h3: " << h3->GetName() << endl;
    return c_new;
}

TCanvas* DrawWithCDF(TH1* h1,
                     TCanvas* c1,
                     std::string CanvasTitle,
                     std::string XaxisTitle,
                     std::string leg_h1,
                     float Xmin,
                     float Xmax,
                     int color) {

    TCanvas* c_new = new TCanvas(CanvasTitle.c_str(), CanvasTitle.c_str(), 800, 600);

    // Define pads
    TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0);
    pad1->SetLeftMargin(0.16);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();

    TPad* pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.315);
    pad2->SetLeftMargin(0.16);pad2->SetBottomMargin(0.33);
    pad2->Draw();

    // Draw upper plot
    pad1->cd();
    c1->DrawClonePad();

    // Clone histogram
    TH1F* h1c = (TH1F*)h1->Clone(TString(h1->GetName()) + "_" + CanvasTitle.c_str() + "_h1");

    // Normalize clone before computing CDF (so CDF goes from 0 to 1)
    if (h1c->Integral() > 0) h1c->Scale(1.0 / h1c->Integral());

    // Build CDF by cumulative sum
    TH1* hCDF1 = h1c->GetCumulative();
    hCDF1->SetName(TString("hCDF1_") + CanvasTitle.c_str());

    // Style CDF1
    hCDF1->SetLineColor(color);
    hCDF1->SetMarkerColor(color);
    hCDF1->SetMarkerStyle(20);

    // Draw CDF in lower pad
    pad2->cd();
    gStyle->SetOptStat(0);
    gPad->SetTickx(0);

    hCDF1->SetTitle("");
    hCDF1->GetYaxis()->SetTitle("CDF");
    hCDF1->GetXaxis()->SetTitle(XaxisTitle.c_str());
    hCDF1->GetYaxis()->SetRangeUser(0, 1);
    hCDF1->GetYaxis()->SetNdivisions(505);
    hCDF1->GetYaxis()->SetTitleFont(43);
    hCDF1->GetXaxis()->SetTitleFont(43);
    hCDF1->GetYaxis()->SetLabelFont(43);
    hCDF1->GetXaxis()->SetLabelFont(43);
    hCDF1->GetYaxis()->SetTitleSize(24);
    hCDF1->GetXaxis()->SetTitleSize(24);
    hCDF1->GetYaxis()->SetLabelSize(20);
    hCDF1->GetXaxis()->SetLabelSize(20);
    hCDF1->GetYaxis()->SetTitleOffset(1.3);
    hCDF1->GetXaxis()->SetTitleOffset(1.0);
    hCDF1->LabelsOption("v", "X");
    hCDF1->GetXaxis()->SetRangeUser(Xmin, Xmax);

    // draw horizontal dashed line at y=0.5
    TLine* line = new TLine(Xmin, 0.5, Xmax, 0.5);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);

    hCDF1->Draw("P");
    line->Draw("same");

    // Legend
    TLegend* leg = new TLegend(0.12, 0.65, 0.45, 0.92);
    leg->AddEntry(hCDF1, leg_h1.c_str(), "pe");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    //leg->Draw();

    c_new->Update();
    cout << "Canvas " << CanvasTitle << " drawn with CDF of: "
         << h1->GetName() << endl;

    return c_new;
}

TCanvas* DrawWithCDF2(TH1* h1,
                      TH1* h2,
                      TCanvas* c1,
                      std::string CanvasTitle,
                      std::string XaxisTitle,
                      std::string leg_h1,
                      std::string leg_h2,
                      float Xmin,
                      float Xmax,
                      int color1,
                      int color2) {

    TCanvas* c_new = new TCanvas(CanvasTitle.c_str(), CanvasTitle.c_str(), 800, 600);

    // Define pads
    TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0);
    pad1->SetLeftMargin(0.16);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();

    TPad* pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.315);
    pad2->SetLeftMargin(0.16);pad2->SetBottomMargin(0.33);
    pad2->Draw();

    // Draw upper plot (the pre-drawn distributions)
    pad1->cd();
    c1->DrawClonePad();

    // Clone histograms
    TH1F* h1c = (TH1F*)h1->Clone(TString(h1->GetName()) + "_" + CanvasTitle.c_str() + "_h1");
    TH1F* h2c = (TH1F*)h2->Clone(TString(h2->GetName()) + "_" + CanvasTitle.c_str() + "_h2");

    // Normalize clones before computing CDF (so CDF goes from 0 to 1)
    if (h1c->Integral() > 0) h1c->Scale(1.0 / h1c->Integral());
    if (h2c->Integral() > 0) h2c->Scale(1.0 / h2c->Integral());

    // Build CDFs by cumulative sum
    TH1* hCDF1 = h1c->GetCumulative();
    hCDF1->SetName(TString("hCDF1_") + CanvasTitle.c_str());
    TH1* hCDF2 = h2c->GetCumulative();
    hCDF2->SetName(TString("hCDF2_") + CanvasTitle.c_str());

    // Style CDF1
    hCDF1->SetLineColor(color1);
    hCDF1->SetMarkerColor(color1);
    hCDF1->SetMarkerStyle(20);

    // Style CDF2
    hCDF2->SetLineColor(color2);
    hCDF2->SetMarkerColor(color2);
    hCDF2->SetMarkerStyle(21);

    // Draw CDFs in lower pad
    pad2->cd();
    gStyle->SetOptStat(0);
    gPad->SetTickx(0);

    hCDF1->SetTitle("");
    hCDF1->GetYaxis()->SetTitle("CDF");
    hCDF1->GetXaxis()->SetTitle(XaxisTitle.c_str());
    hCDF1->GetYaxis()->SetRangeUser(0, 1);
    hCDF1->GetYaxis()->SetNdivisions(505);
    hCDF1->GetYaxis()->SetTitleFont(43);
    hCDF1->GetXaxis()->SetTitleFont(43);
    hCDF1->GetYaxis()->SetLabelFont(43);
    hCDF1->GetXaxis()->SetLabelFont(43);
    hCDF1->GetYaxis()->SetTitleSize(24);
    hCDF1->GetXaxis()->SetTitleSize(24);
    hCDF1->GetYaxis()->SetLabelSize(20);
    hCDF1->GetXaxis()->SetLabelSize(20);
    hCDF1->GetYaxis()->SetTitleOffset(1.3);
    hCDF1->GetXaxis()->SetTitleOffset(1.0);
    hCDF1->LabelsOption("v", "X");
    hCDF1->GetXaxis()->SetRangeUser(Xmin, Xmax);

    // draw horizontal dashed line at y=0.5
    TLine* line = new TLine(Xmin, 0.5, Xmax, 0.5);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);

    hCDF1->Draw("P");
    hCDF2->Draw("P SAME");
    line->Draw("same");

    c_new->Update();
    cout << "Canvas " << CanvasTitle << " drawn with CDFs of: "
         << h1->GetName() << " and " << h2->GetName() << endl;

    return c_new;
}

void ExtractSF (const char *ofiletxt,
                const char *ofiletex,
                const char *htemp_pseudoCaloMET_dataname,
                const char *htemp_orMETtrigger_dataname,
                const char *htemp_pseudoCaloMET_MCname,
                const char *htemp_orMETtrigger_MCname,
                const char *inputfileDATA,
                const char *inputfileMC) {

    TFile *ifileDATA = new TFile(Form("%s", inputfileDATA), "READ");
    TFile *ifileMC   = new TFile(Form("%s", inputfileMC), "READ");

    TH1F *htemp_pseudoCaloMET_data = (TH1F*)ifileDATA->Get(Form("%s", htemp_pseudoCaloMET_dataname));
    TH1F *htemp_orMETtrigger_data = (TH1F*)ifileDATA->Get(Form("%s", htemp_orMETtrigger_dataname));
    TH1F *htemp_pseudoCaloMET_MC   = (TH1F*)ifileMC->Get(Form("%s", htemp_pseudoCaloMET_MCname));
    TH1F *htemp_orMETtrigger_MC   = (TH1F*)ifileMC->Get(Form("%s", htemp_orMETtrigger_MCname));

    ofstream outfile(Form("%s", ofiletxt), ios::out);

    // Construction des edges
    const TAxis* ax = htemp_orMETtrigger_data->GetXaxis();
    const int nBins = ax->GetNbins();

    int nNewBins = 0;
    for (int i = 1; i <= nBins; i++) {
        if (ax->GetBinLowEdge(i) < 400) nNewBins++;
    }
    nNewBins++; // bin [400, +inf]

    double* new_bins = new double[nNewBins + 1];
    int idx = 0;
    for (int i = 1; i <= nBins; i++) {
        if (ax->GetBinLowEdge(i) < 400)
            new_bins[idx++] = ax->GetBinLowEdge(i);
    }
    new_bins[idx++] = 400.0;
    new_bins[idx]   = ax->GetBinUpEdge(nBins);

    // Rebin
    htemp_orMETtrigger_data  = (TH1F*) htemp_orMETtrigger_data ->Rebin(nNewBins, "htemp_orMETtrigger_data",  new_bins);
    htemp_pseudoCaloMET_data = (TH1F*) htemp_pseudoCaloMET_data->Rebin(nNewBins, "htemp_pseudoCaloMET_data", new_bins);
    htemp_orMETtrigger_MC   = (TH1F*) htemp_orMETtrigger_MC   ->Rebin(nNewBins, "htemp_orMETtrigger_MC",    new_bins);
    htemp_pseudoCaloMET_MC  = (TH1F*) htemp_pseudoCaloMET_MC  ->Rebin(nNewBins, "htemp_pseudoCaloMET_MC",   new_bins);
    delete[] new_bins;

    htemp_orMETtrigger_data->Divide(htemp_pseudoCaloMET_data);
    htemp_orMETtrigger_MC->Divide(htemp_pseudoCaloMET_MC);

    htemp_orMETtrigger_data->Divide(htemp_orMETtrigger_MC);
    cout << "Extracted SF for " << ofiletxt << ":\n";
    for (int i = 1; i <= htemp_orMETtrigger_data->GetNbinsX(); i++) {
        double sf, sf_down, sf_up;

        sf      = htemp_orMETtrigger_data->GetBinContent(i);
        sf_down = sf - htemp_orMETtrigger_data->GetBinError(i);
        sf_up   = sf + htemp_orMETtrigger_data->GetBinError(i);

        cout << htemp_orMETtrigger_data->GetBinLowEdge(i) << " " << sf_down << " " << sf << " " << sf_up << "\n";

        outfile << htemp_orMETtrigger_data->GetBinLowEdge(i) << " " << sf_down << " " << sf << " " << sf_up << "\n";
    }
    cout << "\n";
    outfile.close();

    outfile.close();

    // ------------------------------------------------------------------
    // table LaTeX
    // ------------------------------------------------------------------
    std::string textname = std::string(ofiletex);

    std::ofstream tex(textname.c_str(), std::ios::out);

    tex << "\\begin{table}[htbp]\n  \\centering\n"
        << "  \\caption{Trigger scale factors (data/MC) for \\texttt{orMETtrg} "
        << "on the custom binning.}\n"
        << "  \\begin{tabular}{cc}\n    \\hline\n"
        << "    MET bin [GeV] & SF \\\\\n    \\hline\n";

    for (int i = 1; i <= htemp_orMETtrigger_data->GetNbinsX(); i++) {
        double sf  = htemp_orMETtrigger_data->GetBinContent(i);
        double err = htemp_orMETtrigger_data->GetBinError(i);
        double lo  = htemp_orMETtrigger_data->GetBinLowEdge(i);
        double hi  = lo + htemp_orMETtrigger_data->GetBinWidth(i);

        tex << "    $\\left[" << lo << ", " << hi << "\\right]$ & "
            << Form("%.4f $\\pm$ %.4f", sf, err) << " \\\\\n";
    }

    tex << "    \\hline\n  \\end{tabular}\n\\end{table}\n";
    tex.close();


    return;
}

double langaufun(double *x, double *par) {
 
   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.
 
      // Numeric constants
      double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      double mpshift  = -0.22278298;       // Landau maximum location
 
      // Control constants
      double np = 100.0;      // number of convolution steps
      double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
 
      // Variables
      double xx;
      double mpc;
      double fland;
      double sum = 0.0;
      double xlow,xupp;
      double step;
      double i;
 
 
      // MP shift correction
      mpc = par[1] - mpshift * par[0];
 
      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];
 
      step = (xupp-xlow) / np;
 
      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
 
         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }
 
      return (par[2] * step * sum * invsq2pi / par[3]);
}

void PlotNormalized(const char *labelH1, const char *labelH2,
                    const char *inputfile1, const char *inputfile2,
                    const char *histoname1, const char *histoname2,
                    const char *ofilename,
                    const char *xtitle = "x", const char *ytitle = "a.u.",
                    double xmin = 0, double xmax = 1200) {

    gErrorIgnoreLevel = kError;

    TFile *ifile1 = new TFile(inputfile1, "READ");
    TFile *ifile2 = new TFile(inputfile2, "READ");

    if (!ifile1 || ifile1->IsZombie()) {
        std::cerr << "Error: Could not open input file " << inputfile1 << std::endl;
        return;
    }
    if (!ifile2 || ifile2->IsZombie()) {
        std::cerr << "Error: Could not open input file " << inputfile2 << std::endl;
        return;
    }

    TH1F *h1 = (TH1F*)ifile1->Get(Form("%s_%s", labelH1, histoname1));
    TH1F *h2 = (TH1F*)ifile2->Get(Form("%s_%s", labelH2, histoname2));

    if (!h1) {
        std::cerr << "Error: Could not find histo " << Form("%s_%s", labelH1, histoname1) << std::endl;
        return;
    }
    if (!h2) {
        std::cerr << "Error: Could not find histo " << Form("%s_%s", labelH2, histoname2) << std::endl;
        return;
    }

    // detach from files so they survive the TFile closing
    h1->SetDirectory(0);
    h2->SetDirectory(0);

    // normalization to unit area (integral including under/overflow optional)
    if (h1->Integral() > 0) h1->Scale(1.0 / h1->Integral());
    if (h2->Integral() > 0) h2->Scale(1.0 / h2->Integral());

    // styling
    h1->SetLineColor(kRed+1);
    h1->SetMarkerColor(kRed+1);
    h1->SetLineWidth(2);
    h2->SetLineColor(kBlue+1);
    h2->SetMarkerColor(kBlue+1);
    h2->SetLineWidth(2);

    // y range: leave headroom above the tallest bin
    double ymax = std::max(h1->GetMaximum(), h2->GetMaximum());

    TCanvas *c = new TCanvas("c_PlotNormalized", "c_PlotNormalized", 800, 600);
    c->cd();

    h1->GetXaxis()->SetRangeUser(xmin, xmax);
    h1->GetXaxis()->SetTitle(xtitle);
    h1->GetYaxis()->SetTitle(ytitle);
    h1->SetMaximum(1.3 * ymax);
    h1->SetMinimum(1e-8);
    h1->SetStats(0);

    h1->Draw("HIST E");
    h2->Draw("HIST E SAME");

    TLegend *leg = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h1, labelH1, "lep");
    leg->AddEntry(h2, labelH2, "lep");
    leg->Draw();

    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Private work (CMS simulation/data)}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);
    latex1->Draw();

    c->SetLogy();


    c->SaveAs(Form("../outputDisplay/c_PlotNormalized__%s.pdf", ofilename));

    return;
}

TH2D* RebinTH2( const TH2* hIn,
                const std::vector<double>& xBins,
                const std::vector<double>& yBins,
                const TString& newName = "") {
 
    TString hname = newName.IsNull() ?
                    TString(hIn->GetName()) + "_rebinned" :
                    newName;
 
    TH2D* hOut = new TH2D(
        hname,
        hIn->GetTitle(),
        xBins.size() - 1, xBins.data(),
        yBins.size() - 1, yBins.data()
    );
 
    hOut->Sumw2();
 
    for (int ix = 1; ix <= hIn->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= hIn->GetNbinsY(); ++iy) {
 
            double content = hIn->GetBinContent(ix, iy);
            double error   = hIn->GetBinError(ix, iy);
 
            if (content == 0. && error == 0.) continue;
 
            double x = hIn->GetXaxis()->GetBinCenter(ix);
            double y = hIn->GetYaxis()->GetBinCenter(iy);
 
            int binOut = hOut->FindBin(x, y);
 
            double oldContent = hOut->GetBinContent(binOut);
            double oldError   = hOut->GetBinError(binOut);
 
            hOut->SetBinContent(binOut, oldContent + content);
            hOut->SetBinError(
                binOut,
                std::sqrt(oldError * oldError + error * error)
            );
        }
    }
 
    return hOut;
}


// main
void MET_trg_eff(const char* label, const char *ifileName, bool isAOD = false) {

    TFile *ofile;
    if (isAOD) ofile = new TFile("../outputDisplay/MET_trg_eff_AOD.root", "RECREATE");
    else ofile = new TFile("../outputDisplay/MET_trg_eff_miniAOD.root", "RECREATE");

    TFile *ifile = new TFile(Form("%s", ifileName), "READ");

    // input histos
    TH1F *PseudoCaloMET = (TH1F*)ifile->Get(Form("%s_PseudoCaloMET", label));
    TH1F *if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET = (TH1F*)ifile->Get(Form("%s_if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET", label));
    TH1F *if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET = (TH1F*)ifile->Get(Form("%s_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET", label));
    TH1F *if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET = (TH1F*)ifile->Get(Form("%s_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET", label));
    TH1F *if___HLT_MET105_IsoTrk50___PseudoCaloMET = (TH1F*)ifile->Get(Form("%s_if___HLT_MET105_IsoTrk50___PseudoCaloMET", label));
    TH1F *if___orMETtrg___PseudoCaloMET = (TH1F*)ifile->Get(Form("%s_if___orMETtrg___PseudoCaloMET", label));

    TH1F *RecoPFMET = (TH1F*)ifile->Get(Form("%s_RecoPFMET", label));
    TH1F *if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET = (TH1F*)ifile->Get(Form("%s_if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET", label));            
    TH1F *if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET = (TH1F*)ifile->Get(Form("%s_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET", label));
    TH1F *if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET = (TH1F*)ifile->Get(Form("%s_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET", label));
    TH1F *if___HLT_MET105_IsoTrk50___RecoPFMET = (TH1F*)ifile->Get(Form("%s_if___HLT_MET105_IsoTrk50___RecoPFMET", label));
    TH1F *if___orMETtrg___RecoPFMET = (TH1F*)ifile->Get(Form("%s_if___orMETtrg___RecoPFMET", label));
    
    TH1F *PseudoCaloMET__RecoPFMETCut = (TH1F*)ifile->Get(Form("%s_PseudoCaloMET__RecoPFMETCut", label));
    TH1F *if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET__RecoPFMETCut = (TH1F*)ifile->Get(Form("%s_if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET__RecoPFMETCut", label));       
    TH1F *if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET__RecoPFMETCut = (TH1F*)ifile->Get(Form("%s_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET__RecoPFMETCut", label));
    TH1F *if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET__RecoPFMETCut = (TH1F*)ifile->Get(Form("%s_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET__RecoPFMETCut", label));
    TH1F *if___HLT_MET105_IsoTrk50___PseudoCaloMET__RecoPFMETCut = (TH1F*)ifile->Get(Form("%s_if___HLT_MET105_IsoTrk50___PseudoCaloMET__RecoPFMETCut", label));
    TH1F *if___orMETtrg___PseudoCaloMET__RecoPFMETCut = (TH1F*)ifile->Get(Form("%s_if___orMETtrg___PseudoCaloMET__RecoPFMETCut", label));
    
    TH1F *RecoPFMET__PseudoCaloMETCut = (TH1F*)ifile->Get(Form("%s_RecoPFMET__PseudoCaloMETCut", label));
    TH1F *if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET__PseudoCaloMETCut = (TH1F*)ifile->Get(Form("%s_if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET__PseudoCaloMETCut", label));      
    TH1F *if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET__PseudoCaloMETCut = (TH1F*)ifile->Get(Form("%s_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET__PseudoCaloMETCut", label));
    TH1F *if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET__PseudoCaloMETCut = (TH1F*)ifile->Get(Form("%s_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET__PseudoCaloMETCut", label));
    TH1F *if___HLT_MET105_IsoTrk50___RecoPFMET__PseudoCaloMETCut = (TH1F*)ifile->Get(Form("%s_if___HLT_MET105_IsoTrk50___RecoPFMET__PseudoCaloMETCut", label));
    TH1F *if___orMETtrg___RecoPFMET__PseudoCaloMETCut = (TH1F*)ifile->Get(Form("%s_if___orMETtrg___RecoPFMET__PseudoCaloMETCut", label));
        
    TH1F *L1MET = (TH1F*)ifile->Get(Form("%s_L1MET", label));
    TH1F *HLTCaloMET = (TH1F*)ifile->Get(Form("%s_HLTCaloMET", label));
    TH1F *HLTCaloMHT = (TH1F*)ifile->Get(Form("%s_HLTCaloMHT", label));
    TH1F *HLTPFMHT = (TH1F*)ifile->Get(Form("%s_HLTPFMHT", label));
    TH1F *HLTPFMET = (TH1F*)ifile->Get(Form("%s_HLTPFMET", label));


    TH1F *RecoCaloMET;
    TH1F *if___HLT_PFMET120_PFMHT120_IDTight___RecoCaloMET;
    TH1F *if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoCaloMET;
    TH1F *if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoCaloMET;
    TH1F *if___HLT_MET105_IsoTrk50___RecoCaloMET;
    TH1F *if___orMETtrg___RecoCaloMET;
    if (isAOD) {
        RecoCaloMET = (TH1F*)ifile->Get("OnlyMET_RecoCaloMET");
        if___HLT_PFMET120_PFMHT120_IDTight___RecoCaloMET = (TH1F*)ifile->Get("OnlyMET_if___HLT_PFMET120_PFMHT120_IDTight___RecoCaloMET");
        if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoCaloMET = (TH1F*)ifile->Get("OnlyMET_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoCaloMET");
        if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoCaloMET = (TH1F*)ifile->Get("OnlyMET_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoCaloMET");
        if___HLT_MET105_IsoTrk50___RecoCaloMET = (TH1F*)ifile->Get("OnlyMET_if___HLT_MET105_IsoTrk50___RecoCaloMET");
        if___orMETtrg___RecoCaloMET = (TH1F*)ifile->Get("OnlyMET_if___orMETtrg___RecoCaloMET");
    }

    // setup
    TH1F *eff_HLT_PFMET120_PseudoCaloMET = (TH1F*)if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET->Clone("eff_HLT_PFMET120_PseudoCaloMET");
    TH1F *eff_HLT_PFHT500_PseudoCaloMET = (TH1F*)if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET->Clone("eff_HLT_PFHT500_PseudoCaloMET");
    TH1F *eff_HLT_PFMETNoMu120_PseudoCaloMET = (TH1F*)if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET->Clone("eff_HLT_PFMETNoMu120_PseudoCaloMET");
    TH1F *eff_HLT_MET105_PseudoCaloMET = (TH1F*)if___HLT_MET105_IsoTrk50___PseudoCaloMET->Clone("eff_HLT_MET105_PseudoCaloMET");
    TH1F *eff_orMETtrg_PseudoCaloMET = (TH1F*)if___orMETtrg___PseudoCaloMET->Clone("eff_orMETtrg_PseudoCaloMET");
    
    TH1F *eff_HLT_PFMET120_RecoPFMET = (TH1F*)if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET->Clone("eff_HLT_PFMET120_RecoPFMET");
    TH1F *eff_HLT_PFHT500_RecoPFMET = (TH1F*)if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET->Clone("eff_HLT_PFHT500_RecoPFMET");
    TH1F *eff_HLT_PFMETNoMu120_RecoPFMET = (TH1F*)if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET->Clone("eff_HLT_PFMETNoMu120_RecoPFMET");
    TH1F *eff_HLT_MET105_RecoPFMET = (TH1F*)if___HLT_MET105_IsoTrk50___RecoPFMET->Clone("eff_HLT_MET105_RecoPFMET");
    TH1F *eff_orMETtrg_RecoPFMET = (TH1F*)if___orMETtrg___RecoPFMET->Clone("eff_orMETtrg_RecoPFMET");
    
    TH1F *eff_HLT_PFMET120_PseudoCaloMET__RecoPFMETCut = (TH1F*)if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET__RecoPFMETCut->Clone("eff_HLT_PFMET120_PseudoCaloMET__RecoPFMETCut");
    TH1F *eff_HLT_PFHT500_PseudoCaloMET__RecoPFMETCut = (TH1F*)if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET__RecoPFMETCut->Clone("eff_HLT_PFHT500_PseudoCaloMET__RecoPFMETCut");
    TH1F *eff_HLT_PFMETNoMu120_PseudoCaloMET__RecoPFMETCut = (TH1F*)if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET__RecoPFMETCut->Clone("eff_HLT_PFMETNoMu120_PseudoCaloMET__RecoPFMETCut");
    TH1F *eff_HLT_MET105_PseudoCaloMET__RecoPFMETCut = (TH1F*)if___HLT_MET105_IsoTrk50___PseudoCaloMET__RecoPFMETCut->Clone("eff_HLT_MET105_PseudoCaloMET__RecoPFMETCut");
    TH1F *eff_orMETtrg_PseudoCaloMET__RecoPFMETCut = (TH1F*)if___orMETtrg___PseudoCaloMET__RecoPFMETCut->Clone("eff_orMETtrg_PseudoCaloMET__RecoPFMETCut");

    TH1F *eff_HLT_PFMET120_RecoPFMET__PseudoCaloMETCut = (TH1F*)if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET__PseudoCaloMETCut->Clone("eff_HLT_PFMET120_RecoPFMET__PseudoCaloMETCut");
    TH1F *eff_HLT_PFHT500_RecoPFMET__PseudoCaloMETCut = (TH1F*)if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET__PseudoCaloMETCut->Clone("eff_HLT_PFHT500_RecoPFMET__PseudoCaloMETCut");
    TH1F *eff_HLT_PFMETNoMu120_RecoPFMET__PseudoCaloMETCut = (TH1F*)if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET__PseudoCaloMETCut->Clone("eff_HLT_PFMETNoMu120_RecoPFMET__PseudoCaloMETCut");
    TH1F *eff_HLT_MET105_RecoPFMET__PseudoCaloMETCut = (TH1F*)if___HLT_MET105_IsoTrk50___RecoPFMET__PseudoCaloMETCut->Clone("eff_HLT_MET105_RecoPFMET__PseudoCaloMETCut");
    TH1F *eff_orMETtrg_RecoPFMET__PseudoCaloMETCut = (TH1F*)if___orMETtrg___RecoPFMET__PseudoCaloMETCut->Clone("eff_orMETtrg_RecoPFMET__PseudoCaloMETCut");

    TH1F *eff_HLT_PFMET120_RecoCaloMET;
    TH1F *eff_HLT_PFHT500_RecoCaloMET;
    TH1F *eff_HLT_PFMETNoMu120_RecoCaloMET;
    TH1F *eff_HLT_MET105_RecoCaloMET;
    TH1F *eff_orMETtrg_RecoCaloMET;
    if (isAOD) {
        eff_HLT_PFMET120_RecoCaloMET = (TH1F*)if___HLT_PFMET120_PFMHT120_IDTight___RecoCaloMET->Clone("eff_HLT_PFMET120_RecoCaloMET");
        eff_HLT_PFHT500_RecoCaloMET = (TH1F*)if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoCaloMET->Clone("eff_HLT_PFHT500_RecoCaloMET");
        eff_HLT_PFMETNoMu120_RecoCaloMET = (TH1F*)if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoCaloMET->Clone("eff_HLT_PFMETNoMu120_RecoCaloMET");
        eff_HLT_MET105_RecoCaloMET = (TH1F*)if___HLT_MET105_IsoTrk50___RecoCaloMET->Clone("eff_HLT_MET105_RecoCaloMET");
        eff_orMETtrg_RecoCaloMET = (TH1F*)if___orMETtrg___RecoCaloMET->Clone("eff_orMETtrg_RecoCaloMET");
    }

    // efficiency calculation
    eff_HLT_PFMET120_PseudoCaloMET->Divide(PseudoCaloMET);
    eff_HLT_PFHT500_PseudoCaloMET->Divide(PseudoCaloMET);
    eff_HLT_PFMETNoMu120_PseudoCaloMET->Divide(PseudoCaloMET);
    eff_HLT_MET105_PseudoCaloMET->Divide(PseudoCaloMET);
    eff_orMETtrg_PseudoCaloMET->Divide(PseudoCaloMET);
    
    eff_HLT_PFMET120_RecoPFMET->Divide(RecoPFMET);
    eff_HLT_PFHT500_RecoPFMET->Divide(RecoPFMET);
    eff_HLT_PFMETNoMu120_RecoPFMET->Divide(RecoPFMET);
    eff_HLT_MET105_RecoPFMET->Divide(RecoPFMET);
    eff_orMETtrg_RecoPFMET->Divide(RecoPFMET);

    eff_HLT_PFMET120_PseudoCaloMET__RecoPFMETCut->Divide(PseudoCaloMET__RecoPFMETCut);
    eff_HLT_PFHT500_PseudoCaloMET__RecoPFMETCut->Divide(PseudoCaloMET__RecoPFMETCut);
    eff_HLT_PFMETNoMu120_PseudoCaloMET__RecoPFMETCut->Divide(PseudoCaloMET__RecoPFMETCut);
    eff_HLT_MET105_PseudoCaloMET__RecoPFMETCut->Divide(PseudoCaloMET__RecoPFMETCut);
    eff_orMETtrg_PseudoCaloMET__RecoPFMETCut->Divide(PseudoCaloMET__RecoPFMETCut);
    
    eff_HLT_PFMET120_RecoPFMET__PseudoCaloMETCut->Divide(RecoPFMET__PseudoCaloMETCut);
    eff_HLT_PFHT500_RecoPFMET__PseudoCaloMETCut->Divide(RecoPFMET__PseudoCaloMETCut);
    eff_HLT_PFMETNoMu120_RecoPFMET__PseudoCaloMETCut->Divide(RecoPFMET__PseudoCaloMETCut);
    eff_HLT_MET105_RecoPFMET__PseudoCaloMETCut->Divide(RecoPFMET__PseudoCaloMETCut);
    eff_orMETtrg_RecoPFMET__PseudoCaloMETCut->Divide(RecoPFMET__PseudoCaloMETCut);


    if (isAOD) {
        eff_HLT_PFMET120_RecoCaloMET->Divide(RecoCaloMET);
        eff_HLT_PFHT500_RecoCaloMET->Divide(RecoCaloMET);
        eff_HLT_PFMETNoMu120_RecoCaloMET->Divide(RecoCaloMET);
        eff_HLT_MET105_RecoCaloMET->Divide(RecoCaloMET);
        eff_orMETtrg_RecoCaloMET->Divide(RecoCaloMET);
    }


    // drawing

    TCanvas *c_HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET = DrawCanvas(eff_HLT_PFMET120_PseudoCaloMET, "HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET", "Pseudo MET [GeV]", "eff. HLT_PFMET120_PFMHT120_IDTight", "E1", 0, 1200, 0, 1, false);
    TCanvas *c_HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET = DrawCanvas(eff_HLT_PFHT500_PseudoCaloMET, "HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET", "Pseudo MET [GeV]", "eff. HLT_PFHT500_PFMET100_PFMHT100_IDTight", "E1", 0, 1200, 0, 1, false);
    TCanvas *c_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET = DrawCanvas(eff_HLT_PFMETNoMu120_PseudoCaloMET, "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET", "Pseudo MET [GeV]", "eff. HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", "E1", 0, 1200, 0, 1, false);
    TCanvas *c_HLT_MET105_IsoTrk50___PseudoCaloMET = DrawCanvas(eff_HLT_MET105_PseudoCaloMET, "HLT_MET105_IsoTrk50___PseudoCaloMET", "Pseudo MET [GeV]", "eff. HLT_MET105_IsoTrk50", "E1", 0, 1200, 0, 1, false);    
    TCanvas *c_orMETtrg___PseudoCaloMET = DrawCanvas(eff_orMETtrg_PseudoCaloMET, "orMETtrg___PseudoCaloMET", "Pseudo MET [GeV]", "eff. orMETtrg", "E1", 0, 1200, 0, 1, false);

    /*TCanvas *c_HLT_PFMET120_PFMHT120_IDTight___RecoPFMET = DrawCanvas(eff_HLT_PFMET120_RecoPFMET, "HLT_PFMET120_PFMHT120_IDTight___RecoPFMET", "RecoPFMET [GeV]", "eff. HLT_PFMET120_PFMHT120_IDTight", "E1", 0, 1200, 0, 1, false);
    TCanvas *c_HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET = DrawCanvas(eff_HLT_PFHT500_RecoPFMET, "HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET", "RecoPFMET [GeV]", "eff. HLT_PFHT500_PFMET100_PFMHT100_IDTight", "E1", 0, 1200, 0, 1, false);
    TCanvas *c_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET = DrawCanvas(eff_HLT_PFMETNoMu120_RecoPFMET, "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET", "RecoPFMET [GeV]", "eff. HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", "E1", 0, 1200, 0, 1, false);
    TCanvas *c_HLT_MET105_IsoTrk50___RecoPFMET = DrawCanvas(eff_HLT_MET105_RecoPFMET, "HLT_MET105_IsoTrk50___RecoPFMET", "RecoPFMET [GeV]", "eff. HLT_MET105_IsoTrk50", "E1", 0, 1200, 0, 1, false);
    TCanvas *c_orMETtrg___RecoPFMET = DrawCanvas(eff_orMETtrg_RecoPFMET, "orMETtrg___RecoPFMET", "RecoPFMET [GeV]", "eff. orMETtrg", "E1", 0, 1200, 0, 1, false);

    TCanvas *c_HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET__RecoPFMETCut = DrawCanvas(eff_HLT_PFMET120_PseudoCaloMET__RecoPFMETCut, "HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET__RecoPFMETCut", "Pseudo MET [GeV]", "eff. HLT_PFMET120_PFMHT120_IDTight w/RecoPFMET>170 GeV", "E1", 0, 1200, 0, 1, false);
    TCanvas *c_HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET__RecoPFMETCut = DrawCanvas(eff_HLT_PFHT500_PseudoCaloMET__RecoPFMETCut, "HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET__RecoPFMETCut", "Pseudo MET [GeV]", "eff. HLT_PFHT500_PFMET100_PFMHT100_IDTight w/RecoPFMET>170 GeV", "E1", 0, 1200, 0, 1, false);
    TCanvas *c_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET__RecoPFMETCut = DrawCanvas(eff_HLT_PFMETNoMu120_PseudoCaloMET__RecoPFMETCut, "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET__RecoPFMETCut", "Pseudo MET [GeV]", "eff. HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 w/RecoPFMET>170 GeV", "E1", 0, 1200, 0, 1, false);
    TCanvas *c_HLT_MET105_IsoTrk50___PseudoCaloMET__RecoPFMETCut = DrawCanvas(eff_HLT_MET105_PseudoCaloMET__RecoPFMETCut, "HLT_MET105_IsoTrk50___PseudoCaloMET__RecoPFMETCut", "Pseudo MET [GeV]", "eff. HLT_MET105_IsoTrk50 w/RecoPFMET>170 GeV", "E1", 0, 1200, 0, 1, false);
    TCanvas *c_orMETtrg___PseudoCaloMET__RecoPFMETCut = DrawCanvas(eff_orMETtrg_PseudoCaloMET__RecoPFMETCut, "orMETtrg___PseudoCaloMET__RecoPFMETCut", "Pseudo MET [GeV]", "eff. orMETtrgw/RecoPFMET>170 GeV", "E1", 0, 1200, 0, 1, false);

    TCanvas *c_HLT_PFMET120_PFMHT120_IDTight_RecoPFMET__PseudoCaloMETCut = DrawCanvas(eff_HLT_PFMET120_RecoPFMET__PseudoCaloMETCut, "HLT_PFMET120_PFMHT120_IDTight___RecoPFMET__PseudoCaloMETCut", "RecoPFMET [GeV]", "eff. HLT_PFMET120_PFMHT120_IDTight w/Pseudo MET>170 GeV", "E1", 0, 1200, 0, 1, false);
    TCanvas *c_HLT_PFHT500_PFMET100_PFMHT100_IDTight__PseudoCaloMETCut = DrawCanvas(eff_HLT_PFHT500_RecoPFMET__PseudoCaloMETCut, "HLT_PFHT500_PFMET100_PFMHT100_IDTight__PseudoCaloMETCut", "RecoPFMET [GeV]", "eff. HLT_PFHT500_PFMET100_PFMHT100_IDTight w/Pseudo MET>170 GeV", "E1", 0, 1200, 0, 1, false);
    TCanvas *c_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60__PseudoCaloMETCut = DrawCanvas(eff_HLT_PFMETNoMu120_RecoPFMET__PseudoCaloMETCut, "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET__PseudoCaloMETCut", "RecoPFMET [GeV]", "eff. HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 w/Pseudo MET>170 GeV", "E1", 0, 1200, 0, 1, false);
    TCanvas *c_HLT_MET105_IsoTrk50_RecoPFMET__PseudoCaloMETCut = DrawCanvas(eff_HLT_MET105_RecoPFMET__PseudoCaloMETCut, "HLT_MET105_IsoTrk50___RecoPFMET__PseudoCaloMETCut", "RecoPFMET [GeV]", "eff. HLT_MET105_IsoTrk50 w/Pseudo MET>170 GeV", "E1", 0, 1200, 0, 1, false);
    TCanvas *c_orMETtrg_RecoPFMET__PseudoCaloMETCut = DrawCanvas(eff_orMETtrg_RecoPFMET__PseudoCaloMETCut, "orMETtrg_RecoPFMET__PseudoCaloMETCut", "RecoPFMET [GeV]", "eff. orMETtrg w/Pseudo MET>170 GeV", "E1", 0, 1200, 0, 1, false);

    TCanvas *c_PseudoCaloMET = DrawCanvas(PseudoCaloMET, "Pseudo MET", "Pseudo MET [GeV]", "Events / 25 GeV", "HIST", 0, 1000, 0, PseudoCaloMET->GetMaximum()*1.2, true);
    TCanvas *c_RecoPFMET = DrawCanvas(RecoPFMET, "RecoPFMET", "RecoPFMET [GeV]", "Events / 25 GeV", "HIST", 0, 1000, 0, RecoPFMET->GetMaximum()*1.2, true);
    TCanvas *c_L1MET = DrawCanvas(L1MET, "L1MET", "L1MET [GeV]", "Events / 25 GeV", "HIST", 0, 1000, 0, L1MET->GetMaximum()*1.2, true);
    TCanvas *c_HLTCaloMET = DrawCanvas(HLTCaloMET, "HLTCaloMET", "HLTCaloMET [GeV]", "Events / 25 GeV", "HIST", 0, 1000, 0, HLTCaloMET->GetMaximum()*1.2, true);
    TCanvas *c_HLTCaloMHT = DrawCanvas(HLTCaloMHT, "HLTCaloMHT", "HLTCaloMHT [GeV]", "Events / 25 GeV", "HIST", 0, 1000, 0, HLTCaloMHT->GetMaximum()*1.2, true);
    TCanvas *c_HLTPFMHT = DrawCanvas(HLTPFMHT, "HLTPFMHT", "HLTPFMHT [GeV]", "Events / 25 GeV", "HIST", 0, 1000, 0, HLTPFMHT->GetMaximum()*1.2, true);
    TCanvas *c_HLTPFMET = DrawCanvas(HLTPFMET, "HLTPFMET", "HLTPFMET [GeV]", "Events / 25 GeV", "HIST", 0, 1000, 0, HLTPFMET->GetMaximum()*1.2, true);
    */


    // saving
    ofile->cd();
    c_HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET->Write();
    c_HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET->Write();
    c_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET->Write();
    c_HLT_MET105_IsoTrk50___PseudoCaloMET->Write();
    c_orMETtrg___PseudoCaloMET->Write();

    /*c_HLT_PFMET120_PFMHT120_IDTight___RecoPFMET->Write();
    c_HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET->Write();
    c_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET->Write();
    c_HLT_MET105_IsoTrk50___RecoPFMET->Write();
    c_orMETtrg___RecoPFMET->Write();

    c_HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET__RecoPFMETCut->Write();
    c_HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET__RecoPFMETCut->Write();
    c_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET__RecoPFMETCut->Write();
    c_HLT_MET105_IsoTrk50___PseudoCaloMET__RecoPFMETCut->Write();
    c_orMETtrg___PseudoCaloMET__RecoPFMETCut->Write();

    c_HLT_PFMET120_PFMHT120_IDTight_RecoPFMET__PseudoCaloMETCut->Write();
    c_HLT_PFHT500_PFMET100_PFMHT100_IDTight__PseudoCaloMETCut->Write();
    c_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60__PseudoCaloMETCut->Write();
    c_HLT_MET105_IsoTrk50_RecoPFMET__PseudoCaloMETCut->Write();
    c_orMETtrg_RecoPFMET__PseudoCaloMETCut->Write();

    c_PseudoCaloMET->Write();
    c_RecoPFMET->Write();
    c_L1MET->Write();
    c_HLTCaloMET->Write();
    c_HLTCaloMHT->Write();
    c_HLTPFMHT->Write();
    c_HLTPFMET->Write();*/
    
    ofile->Close();


    return;
}


void PFMET_Cut(bool isAOD=false) {

    TFile *ofile = new TFile("../outputDisplay/PFMET_Cut_v2.root", "RECREATE");

    TFile *ifile_AOD;
    if (isAOD) ifile_AOD = new TFile("../output/Gluino2000_Run2_METtrgEff_AOD_V11p15_Eta2p4.root", "READ");
    TFile *ifile_miniAOD = new TFile("../output/Gluino2000_miniAOD_FULL_MET_V11p12_Eta2p4.root", "READ");

    // input histos
    TH1F *RecoPFMET_AOD; TH1F *RecoPFMET_cutPseudoCaloMET_AOD; TH1F *RecoPFMET_cutRecoCaloMET_AOD;
    if (isAOD) {
        RecoPFMET_AOD = (TH1F*)ifile_AOD->Get("RecoPFMET");
        RecoPFMET_cutPseudoCaloMET_AOD = (TH1F*)ifile_AOD->Get("RecoPFMET___wPseudoCaloMETCut");
        RecoPFMET_cutRecoCaloMET_AOD = (TH1F*)ifile_AOD->Get("RecoPFMET___wCaloMETCut");
    }
    

    TH1F *RecoPFMET_miniAOD = (TH1F*)ifile_miniAOD->Get("RecoPFMET");
    TH1F *RecoPFMET_cutPseudoCaloMET_miniAOD = (TH1F*)ifile_miniAOD->Get("RecoPFMET___wPseudoCaloMETCut");

    // setup
    RecoPFMET_miniAOD->SetLineColor(kViolet);
    RecoPFMET_miniAOD->SetMarkerColor(kViolet);
    RecoPFMET_miniAOD->SetMarkerStyle(23);
    RecoPFMET_miniAOD->Rebin(4);
    RecoPFMET_miniAOD->GetXaxis()->SetTitle("RecoPFMET [GeV]");
    RecoPFMET_miniAOD->GetYaxis()->SetTitle("Events");

    if (isAOD) {
        RecoPFMET_cutPseudoCaloMET_AOD->SetLineColor(kRed);
        RecoPFMET_cutPseudoCaloMET_AOD->SetMarkerColor(kRed);
        RecoPFMET_cutPseudoCaloMET_AOD->SetMarkerStyle(20);
        RecoPFMET_cutPseudoCaloMET_AOD->GetXaxis()->SetTitle("RecoPFMET [GeV]");
        RecoPFMET_cutPseudoCaloMET_AOD->GetYaxis()->SetTitle("Events");
        RecoPFMET_cutPseudoCaloMET_AOD->Rebin(4);

        RecoPFMET_AOD->SetLineColor(kBlack);
        RecoPFMET_AOD->SetMarkerColor(kBlack);
        RecoPFMET_AOD->SetMarkerStyle(43);
        RecoPFMET_AOD->Rebin(4);
        RecoPFMET_AOD->GetXaxis()->SetTitle("RecoPFMET [GeV]");
        RecoPFMET_AOD->GetYaxis()->SetTitle("Events");

        RecoPFMET_cutRecoCaloMET_AOD->SetLineColor(kGreen+2);
        RecoPFMET_cutRecoCaloMET_AOD->SetMarkerColor(kGreen+2);
        RecoPFMET_cutRecoCaloMET_AOD->SetMarkerStyle(22);
        RecoPFMET_cutRecoCaloMET_AOD->Rebin(4);
    }

    RecoPFMET_cutPseudoCaloMET_miniAOD->SetLineColor(kBlue);
    RecoPFMET_cutPseudoCaloMET_miniAOD->SetMarkerColor(kBlue);
    RecoPFMET_cutPseudoCaloMET_miniAOD->SetMarkerStyle(21);
    RecoPFMET_cutPseudoCaloMET_miniAOD->Rebin(4);

    


    // drawing
    TCanvas *c_PFMET_AODcut_ratio; TCanvas *c_PFMET_nocut_ratio;
    if (isAOD) {
        cout << "c_PFMET_nocut" << endl;
        TCanvas *c_PFMET_nocut = new TCanvas("c_PFMET_nocut","c_PFMET_nocut",800,800);
        c_PFMET_nocut->cd();
        RecoPFMET_AOD->Draw("E1");
        RecoPFMET_miniAOD->Draw("E1 same");
        TLegend *leg_PFMET_nocut = new TLegend(0.6,0.7,0.88,0.88);
        leg_PFMET_nocut->AddEntry(RecoPFMET_AOD, "AOD", "ep");
        leg_PFMET_nocut->AddEntry(RecoPFMET_miniAOD, "miniAOD", "ep");
        leg_PFMET_nocut->Draw("same");
        leg_PFMET_nocut->SetBorderSize(0);
        leg_PFMET_nocut->SetFillStyle(0);
        //c_PFMET_nocut_ratio = DrawWithRatio(RecoPFMET_AOD, RecoPFMET_miniAOD, c_PFMET_nocut, "RecoPFMET AOD vs miniAOD", "AOD/miniAOD", "E1 same", "E1 same", 0, 2000);

        cout << "c_PFMET_AODcut" << endl;
        TCanvas *c_PFMET_AODcut = new TCanvas("c_PFMET_AODcut","c_PFMET_AODcut",800,800);
        c_PFMET_AODcut->cd();
        RecoPFMET_AOD->Draw("E1");
        RecoPFMET_cutPseudoCaloMET_AOD->Draw("E1 same");
        TLegend *leg_PFMET_AODcut = new TLegend(0.6,0.7,0.88,0.88);
        leg_PFMET_AODcut->AddEntry(RecoPFMET_AOD, "AOD wo PseudoCaloMET cut", "ep");
        leg_PFMET_AODcut->AddEntry(RecoPFMET_cutPseudoCaloMET_AOD, "AOD with PseudoCaloMET cut", "ep");
        leg_PFMET_AODcut->Draw("same");
        leg_PFMET_AODcut->SetBorderSize(0);
        leg_PFMET_AODcut->SetFillStyle(0);
        //c_PFMET_AODcut_ratio = DrawWithRatio(RecoPFMET_AOD, RecoPFMET_cutPseudoCaloMET_AOD, c_PFMET_AODcut, "RecoPFMET AOD cut", "wo cut/cut", "E1 same", "E1 same", 0, 2000);
    }

    cout << "c_PFMET_miniAODcut" << endl;
    TCanvas *c_PFMET_miniAODcut = new TCanvas("c_PFMET_miniAODcut","c_PFMET_miniAODcut",800,800);
    c_PFMET_miniAODcut->cd();
    RecoPFMET_miniAOD->Draw("E1");
    RecoPFMET_cutPseudoCaloMET_miniAOD->Draw("E1 same");
    TLegend *leg_PFMET_miniAODcut = new TLegend(0.6,0.7,0.88,0.88);
    leg_PFMET_miniAODcut->AddEntry(RecoPFMET_miniAOD, "miniAOD wo PseudoCaloMET cut", "ep");
    leg_PFMET_miniAODcut->AddEntry(RecoPFMET_cutPseudoCaloMET_miniAOD, "miniAOD with PseudoCaloMET cut", "ep");
    leg_PFMET_miniAODcut->Draw("same");
    leg_PFMET_miniAODcut->SetBorderSize(0);
    leg_PFMET_miniAODcut->SetFillStyle(0);
    //TCanvas *c_PFMET_miniAODcut_ratio = DrawWithRatio(RecoPFMET_miniAOD, RecoPFMET_cutPseudoCaloMET_miniAOD, c_PFMET_miniAODcut, "RecoPFMET miniAOD cut", "wo cut/cut", "E1 same", "E1 same", 0, 2000);


    TCanvas *c_PFMET_CaloMETcut_ratio; TCanvas *c_PFMET_ratio;
    if (isAOD) {
        cout << "c_PFMET_CaloMETcut" << endl;
        TCanvas *c_PFMET_CaloMETcut = new TCanvas("c_PFMET_CaloMETcut","c_PFMET_CaloMETcut",800,800);
        c_PFMET_CaloMETcut->cd();
        RecoPFMET_cutPseudoCaloMET_AOD->Draw("E1");
        RecoPFMET_cutRecoCaloMET_AOD->Draw("E1 same");
        TLegend *leg_PFMET_CaloMETcut = new TLegend(0.6,0.7,0.88,0.88);
        leg_PFMET_CaloMETcut->AddEntry(RecoPFMET_cutPseudoCaloMET_AOD, "AOD w PseudoCaloMET cut", "ep");
        leg_PFMET_CaloMETcut->AddEntry(RecoPFMET_cutRecoCaloMET_AOD, "AOD w RecoCaloMET cut", "ep");
        leg_PFMET_CaloMETcut->Draw("same");
        leg_PFMET_CaloMETcut->SetBorderSize(0);
        leg_PFMET_CaloMETcut->SetFillStyle(0);
        //c_PFMET_CaloMETcut_ratio = DrawWithRatio(RecoPFMET_cutPseudoCaloMET_AOD, RecoPFMET_cutRecoCaloMET_AOD, c_PFMET_CaloMETcut, "RecoPFMET AOD", "PseudoCaloMET cut/RecoCaloMET cut", "E1 same", "E1 same", 0, 2000);
    
            cout << "c_PFMET_onlycut" << endl;
        TCanvas *c_PFMET_onlycut = new TCanvas("c_PFMET_onlycut","c_PFMET_onlycut",800,800);
        c_PFMET_onlycut->cd();
        RecoPFMET_cutPseudoCaloMET_AOD->Draw("E1");
        RecoPFMET_cutPseudoCaloMET_miniAOD->Draw("E1 same");
        TLegend *leg_PFMET = new TLegend(0.6,0.7,0.88,0.88);
        leg_PFMET->AddEntry(RecoPFMET_cutPseudoCaloMET_AOD, "AOD", "ep");
        leg_PFMET->AddEntry(RecoPFMET_cutPseudoCaloMET_miniAOD, "miniAOD", "ep");
        leg_PFMET->Draw("same");
        leg_PFMET->SetBorderSize(0);
        leg_PFMET->SetFillStyle(0);
        //c_PFMET_ratio = DrawWithRatio(RecoPFMET_cutPseudoCaloMET_AOD, RecoPFMET_cutPseudoCaloMET_miniAOD, c_PFMET_onlycut, "RecoPFMET with PseudoCaloMET cut", "AOD/miniAOD", "E1 same", "E1 same", 0, 2000);
    }
    // saving
    cout << "saving histos..." << endl;
    ofile->cd();
    //c_PFMET_miniAODcut_ratio->Write();
    if (isAOD) {
        c_PFMET_CaloMETcut_ratio->Write();
        c_PFMET_AODcut_ratio->Write();
        c_PFMET_ratio->Write();
        c_PFMET_nocut_ratio->Write();
    }
    ofile->Close();

    return;

}


void TrigEff_AODvsMiniAOD() {

    TFile *ofile = new TFile("../outputDisplay/TrigEff_AODvsMiniAOD.root", "RECREATE");

    //TFile *ifile_AOD = new TFile("../output/Gluino2000_Run2_METtrgEff_AOD_V11p15_Eta2p4.root", "READ");
    //TFile *ifile_miniAOD = new TFile("../output/Gluino2000_Run2_METtrgEff_V11p15_Eta2p4.root", "READ");

    TFile *ifile_AOD = new TFile("../output/Gluino2000_Run2_MET_AOD_V11p16_Eta2p4.root", "READ");
    TFile *ifile_miniAOD = new TFile("../output/Gluino2000_Run2_MET_V11p16_Eta2p4.root", "READ");


    // input histos
    TH1F *RecoPFMET_AOD = (TH1F*)ifile_AOD->Get("OnlyMET_RecoPFMET");
    TH1F *RecoPFMET_miniAOD = (TH1F*)ifile_miniAOD->Get("METanalysis_RecoPFMET");
    TH1F *if___orMETtrg___RecoPFMET_AOD = (TH1F*)ifile_AOD->Get("OnlyMET_if___orMETtrg___RecoPFMET");
    TH1F *if___orMETtrg___RecoPFMET_miniAOD = (TH1F*)ifile_miniAOD->Get("METanalysis_if___orMETtrg___RecoPFMET");

    TH1F *PseudoCaloMET_AOD = (TH1F*)ifile_AOD->Get("OnlyMET_CaloJets");
    TH1F *PseudoCaloMET_miniAOD = (TH1F*)ifile_miniAOD->Get("METanalysis_CaloJets");
    TH1F *if___orMETtrg___PseudoCaloMET_AOD = (TH1F*)ifile_AOD->Get("OnlyMET_if___orMETtrg___CaloJets");
    TH1F *if___orMETtrg___PseudoCaloMET_miniAOD = (TH1F*)ifile_miniAOD->Get("METanalysis_if___orMETtrg___CaloJets");

    
    // setup
    TH1F *eff_orMETtrg_RecoPFMET_AOD = (TH1F*)if___orMETtrg___RecoPFMET_AOD->Clone("eff_orMETtrg_RecoPFMET_AOD");
    TH1F *eff_orMETtrg_RecoPFMET_miniAOD = (TH1F*)if___orMETtrg___RecoPFMET_miniAOD->Clone("eff_orMETtrg_RecoPFMET_miniAOD");
    
    TH1F *eff_orMETtrg_PseudoCaloMET_AOD = (TH1F*)if___orMETtrg___PseudoCaloMET_AOD->Clone("eff_orMETtrg_PseudoCaloMET_AOD");
    TH1F *eff_orMETtrg_PseudoCaloMET_miniAOD = (TH1F*)if___orMETtrg___PseudoCaloMET_miniAOD->Clone("eff_orMETtrg_PseudoCaloMET_miniAOD");

    eff_orMETtrg_RecoPFMET_AOD->Divide(RecoPFMET_AOD);
    eff_orMETtrg_RecoPFMET_miniAOD->Divide(RecoPFMET_miniAOD);
    eff_orMETtrg_PseudoCaloMET_AOD->Divide(PseudoCaloMET_AOD);
    eff_orMETtrg_PseudoCaloMET_miniAOD->Divide(PseudoCaloMET_miniAOD);


    // styling
    eff_orMETtrg_RecoPFMET_AOD->SetLineColor(kRed);
    eff_orMETtrg_RecoPFMET_AOD->SetMarkerColor(kRed);
    eff_orMETtrg_RecoPFMET_AOD->SetMarkerStyle(21);
    eff_orMETtrg_RecoPFMET_AOD->GetXaxis()->SetTitle("RecoPFMET [GeV]");
    eff_orMETtrg_RecoPFMET_AOD->GetYaxis()->SetTitle("eff. orMETtrg");
    eff_orMETtrg_RecoPFMET_AOD->GetXaxis()->SetRangeUser(0, 1500);
    eff_orMETtrg_RecoPFMET_AOD->GetYaxis()->SetRangeUser(0, 1);

    eff_orMETtrg_RecoPFMET_miniAOD->SetLineColor(kBlue);
    eff_orMETtrg_RecoPFMET_miniAOD->SetMarkerColor(kBlue);
    eff_orMETtrg_RecoPFMET_miniAOD->SetMarkerStyle(22);
    eff_orMETtrg_RecoPFMET_miniAOD->GetXaxis()->SetTitle("RecoPFMET [GeV]");
    eff_orMETtrg_RecoPFMET_miniAOD->GetYaxis()->SetTitle("eff. orMETtrg");
    eff_orMETtrg_RecoPFMET_miniAOD->GetXaxis()->SetRangeUser(0, 1500);
    eff_orMETtrg_RecoPFMET_miniAOD->GetYaxis()->SetRangeUser(0, 1);

    eff_orMETtrg_PseudoCaloMET_AOD->SetLineColor(kRed);
    eff_orMETtrg_PseudoCaloMET_AOD->SetMarkerColor(kRed);
    eff_orMETtrg_PseudoCaloMET_AOD->SetMarkerStyle(21);
    eff_orMETtrg_PseudoCaloMET_AOD->GetXaxis()->SetTitle("PseudoMET [GeV]");
    eff_orMETtrg_PseudoCaloMET_AOD->GetYaxis()->SetTitle("eff. orMETtrg");
    eff_orMETtrg_PseudoCaloMET_AOD->GetXaxis()->SetRangeUser(0, 1500);
    eff_orMETtrg_PseudoCaloMET_AOD->GetYaxis()->SetRangeUser(0, 1);

    eff_orMETtrg_PseudoCaloMET_miniAOD->SetLineColor(kBlue);
    eff_orMETtrg_PseudoCaloMET_miniAOD->SetMarkerColor(kBlue);
    eff_orMETtrg_PseudoCaloMET_miniAOD->SetMarkerStyle(22);
    eff_orMETtrg_PseudoCaloMET_miniAOD->GetXaxis()->SetTitle("PseudoMET [GeV]");
    eff_orMETtrg_PseudoCaloMET_miniAOD->GetYaxis()->SetTitle("eff. orMETtrg");
    eff_orMETtrg_PseudoCaloMET_miniAOD->GetXaxis()->SetRangeUser(0, 1500);
    eff_orMETtrg_PseudoCaloMET_miniAOD->GetYaxis()->SetRangeUser(0, 1);


    // drawing
    TCanvas *c_PFMET = new TCanvas("c_PFMET","c_PFMET",800,800);
    c_PFMET->cd();
    eff_orMETtrg_RecoPFMET_miniAOD->Draw("E1");
    eff_orMETtrg_RecoPFMET_AOD->Draw("E1 same");
    TLegend *leg_PFMET = new TLegend(0.6,0.7,0.88,0.88);
    leg_PFMET->AddEntry(eff_orMETtrg_RecoPFMET_AOD, "AOD", "ep");
    leg_PFMET->AddEntry(eff_orMETtrg_RecoPFMET_miniAOD, "miniAOD", "ep");
    leg_PFMET->Draw("same");
    leg_PFMET->SetBorderSize(0);
    leg_PFMET->SetFillStyle(0);
    //TCanvas *c_PFMET_ratio = DrawWithRatio(eff_orMETtrg_RecoPFMET_miniAOD, eff_orMETtrg_RecoPFMET_AOD, c_PFMET, "RecoPFMET", "miniAOD/AOD", "E1 same", "E1 same", 0, 1200, 0, 1, false);

    TCanvas *c_PseudoCaloMET = new TCanvas("c_PseudoCaloMET","c_PseudoCaloMET",800,800);
    c_PseudoCaloMET->cd();
    eff_orMETtrg_PseudoCaloMET_miniAOD->Draw("E1");
    eff_orMETtrg_PseudoCaloMET_AOD->Draw("E1 same");
    TLegend *leg_PseudoCaloMET = new TLegend(0.6,0.7,0.88,0.88);
    leg_PseudoCaloMET->AddEntry(eff_orMETtrg_PseudoCaloMET_AOD, "AOD", "ep");
    leg_PseudoCaloMET->AddEntry(eff_orMETtrg_PseudoCaloMET_miniAOD, "miniAOD", "ep");
    leg_PseudoCaloMET->Draw("same");
    leg_PseudoCaloMET->SetBorderSize(0);
    leg_PseudoCaloMET->SetFillStyle(0);
    //TCanvas *c_PseudoCaloMET_ratio = DrawWithRatio(eff_orMETtrg_PseudoCaloMET_miniAOD, eff_orMETtrg_PseudoCaloMET_AOD, c_PseudoCaloMET, "PseudoCaloMET", "miniAOD/AOD", "E1 same", "E1 same", 0, 1200, 0, 1, false);


    // saving
    ofile->cd();
    c_PFMET->Write();
    //c_PFMET_ratio->Write();
    c_PseudoCaloMET->Write();
    //c_PseudoCaloMET_ratio->Write();
    ofile->Close();

    return;
}


void Cutflows(std::string QCD,
              std::string TTbar,
              std::string TTbarSemiLep,
              std::string Wjets,
              std::string JetMETdata,
              std::string Gluino2000) {

    TFile *ofile = new TFile(OutPath("EventCutflow.root"), "RECREATE");

    TFile *ifile_JetMETdata = new TFile(JetMETdata.c_str(), "READ");
    TFile *ifile_Gluino2000 = new TFile(Gluino2000.c_str(), "READ");
    TH1D *Cutflow_JetMETdata = (TH1D*)ifile_JetMETdata->Get("EventCutflow");
    TH1D *Cutflow_Gluino2000 = (TH1D*)ifile_Gluino2000->Get("EventCutflow");

    // backgrounds in stacking order; an empty path leaves the sample out
    struct Bkg { std::string path; const char *key; const char *label; int color; TH1D *h; };
    std::vector<Bkg> bkgs;
    for (Bkg b : std::vector<Bkg>{{QCD, "QCD", "QCD (#mu enriched)", kGreen, nullptr},
                                  {TTbar, "TTbar", "t#bar{t}#rightarrow2l2#nu", kRed, nullptr},
                                  {TTbarSemiLep, "TTbarSemiLep", "t#bar{t}#rightarrowl#nu", kRed+2, nullptr},
                                  {Wjets, "Wjets", "W(#rightarrow#mu#nu)+jets", kBlue-7, nullptr}}) {
        TFile *f = OpenOptional(b.path.c_str());
        if (!f) { std::cout << "Cutflows: " << b.key << " left out" << std::endl; continue; }
        b.h = (TH1D*)f->Get("EventCutflow");
        if (!b.h) Abort(Form("no EventCutflow in %s", b.path.c_str()));
        bkgs.push_back(b);
    }
    if (bkgs.empty()) Abort("Cutflows: give at least one background");

    THStack *hs = new THStack("hs", "");
    for (auto &b : bkgs) hs->Add(b.h);

    auto PrintCutflowPercent = [](TH1D* h, const std::string& name, int cut = 1)
    {
        double N0 = h->GetBinContent(cut);

        std::cout << "\n=== " << name << " ===\n";
        std::cout << "Cut\t\tEvents\t\tRemaining (%)\n";

        for (int i = 2; i <= h->GetNbinsX(); ++i) { // remove 'All' bin
            double Ni = h->GetBinContent(i);
            double eff = (N0 > 0) ? 100. * Ni / N0 : 0.;
            std::cout << i << "\t\t"
                    << Ni << "\t\t"
                    << std::fixed << std::setprecision(3)
                    << eff << " %\n";
        }
    };

    // --- Style   
    int colorSignal = kOrange+8;
    int colorJetMET = kBlack;
    gStyle->SetOptStat(0);

    for (auto &b : bkgs) {
        b.h->SetFillColorAlpha(b.color, 0.5);
        b.h->SetLineColor(b.color);
    }
    Cutflow_Gluino2000->SetLineColor(colorSignal);
    Cutflow_Gluino2000->SetMarkerColor(colorSignal);
    Cutflow_Gluino2000->SetMarkerStyle(22);
    Cutflow_JetMETdata->SetLineColor(colorJetMET);
    Cutflow_JetMETdata->SetMarkerColor(colorJetMET);
    Cutflow_JetMETdata->SetMarkerStyle(20);

    TLatex *tex = new TLatex(0.68, 0.91, Label("lumi", "109 fb^{-1} (13.6 TeV)"));
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);

    TLatex *latex1 = new TLatex(0.16, 0.91, "#it{Private work (CMS simulation/data)}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);

    Cutflow_Gluino2000->Scale(100);


    // --- Canvas + pads
    TCanvas* c = new TCanvas("Cutflow_BKG_data_signal", "", 800, 600);
    c->SetBottomMargin(0.28);
    c->SetLeftMargin(0.16);
    c->cd();

    hs->SetMaximum(1e9);
    hs->SetMinimum(5e2);

    hs->Draw("HIST");
    hs->GetHistogram()->GetYaxis()->SetTitle("Number of events");
    Cutflow_Gluino2000->Draw("E1 same");
    Cutflow_JetMETdata->Draw("E1 same");
    c->SetLogy();

    // --- Légende
    TLegend* leg = new TLegend(0.50, 0.72, 0.89, 0.89);
    leg->SetBorderSize(0);
    leg->SetNColumns(2);
    leg->AddEntry(Cutflow_JetMETdata, Label("JetMETdata", "MET data"),"lep");
    for (const char *key : {"Wjets", "TTbar", "TTbarSemiLep", "QCD"})      // legend order
        for (auto &b : bkgs)
            if (!strcmp(b.key, key)) leg->AddEntry(b.h, Label(b.key, b.label), "f");
    leg->AddEntry(Cutflow_Gluino2000, Form("%s x100", Label("Gluino2000", "#tilde{g} (m=2000 GeV)")), "ep");
    leg->Draw();

    // --- Label CMS
    latex1->Draw();
    tex->Draw();

    // --- Labels X
    std::vector<TString> xlabel = {"All","or MET trigger","MET filters", "PUppi MET > 150", "p_{T}^{PF and track}>50","|#eta|<2.4","N_{pixel hits}#geq2","N_{dEdx hits}#geq10", "f_{valid hits}>0.8"
                ,"High Purity","#chi^{2}/N_{dof}<5","|d_{z}|<0.1","|d_{xy}|<0.02","I^{rel}_{PF}<0.02","I^{trk}_{dr03}<15",
                "E/p<0.3","#sigma_{p_{T}}/p_{T}^{2}<0.0008","#sigma_{p_{T}}/p_{T}<1","F_{pixel}>0.3","I_{h}>C"};

    for (int i = 1; i <= bkgs[0].h->GetNbinsX(); i++) {
        if (i-1 < xlabel.size()) {
            for (auto &b : bkgs) b.h->GetXaxis()->SetBinLabel(i, xlabel[i-1]);
            hs->GetXaxis()->SetBinLabel(i, xlabel[i-1]);
        }
    }

    gPad->SetTickx(0);
    hs->GetYaxis()->SetTitleSize(0.06);
    hs->GetYaxis()->SetTitleOffset(0.9);
    hs->GetYaxis()->SetLabelSize(0.05);
    hs->GetXaxis()->LabelsOption("v");
    hs->GetXaxis()->SetRangeUser(1, 20);
    hs->GetXaxis()->SetLabelSize(0.05);

    gPad->SetTicky(1);
    gPad->SetTickx(1);

    c->Update();


    // Cout values:
    PrintCutflowPercent(Cutflow_JetMETdata, "JetMET data", 1);
    for (auto &b : bkgs) PrintCutflowPercent(b.h, b.key, 1);
    PrintCutflowPercent(Cutflow_Gluino2000, "Signal m=2000", 1);



    // --- Saving
    ofile->cd();
    c->Write();
    c->SaveAs(OutPath("Nm1plots/EventCutflow.pdf"));
    ofile->Close();
}



void MET_trg_eff(const char *labelData, 
                 const char *labelMC, 
                 const char *ofilename, 
                 const char *inputfileDATA, 
                 const char *inputfileMC) {

    cout << "dataset: " << inputfileDATA << " and " << inputfileMC << endl;

    gErrorIgnoreLevel = kError;

    TFile *ofile = new TFile(Form("../outputDisplay/%s.root", ofilename), "RECREATE");

    TFile *ifileDATA = new TFile(Form("%s", inputfileDATA), "READ");
    TFile *ifileMC = new TFile(Form("%s", inputfileMC), "READ");


    // MC 
    TH1F *PseudoCaloMET_MC = (TH1F*)ifileMC->Get(Form("%s_PseudoCaloMET", labelMC));
    TH1F *if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET_MC = (TH1F*)ifileMC->Get(Form("%s_if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET", labelMC));
    TH1F *if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET_MC = (TH1F*)ifileMC->Get(Form("%s_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET", labelMC));
    TH1F *if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET_MC = (TH1F*)ifileMC->Get(Form("%s_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET", labelMC));
    TH1F *if___HLT_MET105_IsoTrk50___PseudoCaloMET_MC = (TH1F*)ifileMC->Get(Form("%s_if___HLT_MET105_IsoTrk50___PseudoCaloMET", labelMC));
    TH1F *if___orMETtrg___PseudoCaloMET_MC = (TH1F*)ifileMC->Get(Form("%s_if___orMETtrg___PseudoCaloMET", labelMC));

    TH1F *RecoPFMET_MC = (TH1F*)ifileMC->Get(Form("%s_RecoPFMET", labelMC));
    TH1F *if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET_MC = (TH1F*)ifileMC->Get(Form("%s_if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET", labelMC));            
    TH1F *if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET_MC = (TH1F*)ifileMC->Get(Form("%s_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET", labelMC));
    TH1F *if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET_MC = (TH1F*)ifileMC->Get(Form("%s_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET", labelMC));
    TH1F *if___HLT_MET105_IsoTrk50___RecoPFMET_MC = (TH1F*)ifileMC->Get(Form("%s_if___HLT_MET105_IsoTrk50___RecoPFMET", labelMC));
    TH1F *if___orMETtrg___RecoPFMET_MC = (TH1F*)ifileMC->Get(Form("%s_if___orMETtrg___RecoPFMET", labelMC));
    
    TH1F *PseudoCaloMET__RecoPFMETCut_MC = (TH1F*)ifileMC->Get(Form("%s_PseudoCaloMET__RecoPFMETCut", labelMC));
    TH1F *if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET__RecoPFMETCut_MC = (TH1F*)ifileMC->Get(Form("%s_if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET__RecoPFMETCut", labelMC));       
    TH1F *if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET__RecoPFMETCut_MC = (TH1F*)ifileMC->Get(Form("%s_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET__RecoPFMETCut", labelMC));
    TH1F *if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET__RecoPFMETCut_MC = (TH1F*)ifileMC->Get(Form("%s_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET__RecoPFMETCut", labelMC));
    TH1F *if___HLT_MET105_IsoTrk50___PseudoCaloMET__RecoPFMETCut_MC = (TH1F*)ifileMC->Get(Form("%s_if___HLT_MET105_IsoTrk50___PseudoCaloMET__RecoPFMETCut", labelMC));
    TH1F *if___orMETtrg___PseudoCaloMET__RecoPFMETCut_MC = (TH1F*)ifileMC->Get(Form("%s_if___orMETtrg___PseudoCaloMET__RecoPFMETCut", labelMC));
    
    TH1F *RecoPFMET__PseudoCaloMETCut_MC = (TH1F*)ifileMC->Get(Form("%s_RecoPFMET__PseudoCaloMETCut", labelMC));
    TH1F *if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET__PseudoCaloMETCut_MC = (TH1F*)ifileMC->Get(Form("%s_if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET__PseudoCaloMETCut", labelMC));      
    TH1F *if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET__PseudoCaloMETCut_MC = (TH1F*)ifileMC->Get(Form("%s_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET__PseudoCaloMETCut", labelMC));
    TH1F *if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET__PseudoCaloMETCut_MC = (TH1F*)ifileMC->Get(Form("%s_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET__PseudoCaloMETCut", labelMC));
    TH1F *if___HLT_MET105_IsoTrk50___RecoPFMET__PseudoCaloMETCut_MC = (TH1F*)ifileMC->Get(Form("%s_if___HLT_MET105_IsoTrk50___RecoPFMET__PseudoCaloMETCut", labelMC));
    TH1F *if___orMETtrg___RecoPFMET__PseudoCaloMETCut_MC = (TH1F*)ifileMC->Get(Form("%s_if___orMETtrg___RecoPFMET__PseudoCaloMETCut", labelMC));


    // DATA
    TH1F *PseudoCaloMET_DATA = (TH1F*)ifileDATA->Get(Form("%s_PseudoCaloMET", labelData));
    TH1F *if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET_DATA = (TH1F*)ifileDATA->Get(Form("%s_if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET", labelData));
    TH1F *if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET_DATA = (TH1F*)ifileDATA->Get(Form("%s_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET", labelData));
    TH1F *if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET_DATA = (TH1F*)ifileDATA->Get(Form("%s_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET", labelData));
    TH1F *if___HLT_MET105_IsoTrk50___PseudoCaloMET_DATA = (TH1F*)ifileDATA->Get(Form("%s_if___HLT_MET105_IsoTrk50___PseudoCaloMET", labelData));
    TH1F *if___orMETtrg___PseudoCaloMET_DATA = (TH1F*)ifileDATA->Get(Form("%s_if___orMETtrg___PseudoCaloMET", labelData));

    TH1F *RecoPFMET_DATA = (TH1F*)ifileDATA->Get(Form("%s_RecoPFMET", labelData));
    TH1F *if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET_DATA = (TH1F*)ifileDATA->Get(Form("%s_if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET", labelData));            
    TH1F *if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET_DATA = (TH1F*)ifileDATA->Get(Form("%s_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET", labelData));
    TH1F *if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET_DATA = (TH1F*)ifileDATA->Get(Form("%s_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET", labelData));
    TH1F *if___HLT_MET105_IsoTrk50___RecoPFMET_DATA = (TH1F*)ifileDATA->Get(Form("%s_if___HLT_MET105_IsoTrk50___RecoPFMET", labelData));
    TH1F *if___orMETtrg___RecoPFMET_DATA = (TH1F*)ifileDATA->Get(Form("%s_if___orMETtrg___RecoPFMET", labelData));
    
    TH1F *PseudoCaloMET__RecoPFMETCut_DATA = (TH1F*)ifileDATA->Get(Form("%s_PseudoCaloMET__RecoPFMETCut", labelData));
    TH1F *if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET__RecoPFMETCut_DATA = (TH1F*)ifileDATA->Get(Form("%s_if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET__RecoPFMETCut", labelData));       
    TH1F *if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET__RecoPFMETCut_DATA = (TH1F*)ifileDATA->Get(Form("%s_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET__RecoPFMETCut", labelData));
    TH1F *if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET__RecoPFMETCut_DATA = (TH1F*)ifileDATA->Get(Form("%s_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET__RecoPFMETCut", labelData));
    TH1F *if___HLT_MET105_IsoTrk50___PseudoCaloMET__RecoPFMETCut_DATA = (TH1F*)ifileDATA->Get(Form("%s_if___HLT_MET105_IsoTrk50___PseudoCaloMET__RecoPFMETCut", labelData));
    TH1F *if___orMETtrg___PseudoCaloMET__RecoPFMETCut_DATA = (TH1F*)ifileDATA->Get(Form("%s_if___orMETtrg___PseudoCaloMET__RecoPFMETCut", labelData));
    
    TH1F *RecoPFMET__PseudoCaloMETCut_DATA = (TH1F*)ifileDATA->Get(Form("%s_RecoPFMET__PseudoCaloMETCut", labelData));
    TH1F *if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET__PseudoCaloMETCut_DATA = (TH1F*)ifileDATA->Get(Form("%s_if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET__PseudoCaloMETCut", labelData));      
    TH1F *if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET__PseudoCaloMETCut_DATA = (TH1F*)ifileDATA->Get(Form("%s_if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET__PseudoCaloMETCut", labelData));
    TH1F *if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET__PseudoCaloMETCut_DATA = (TH1F*)ifileDATA->Get(Form("%s_if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET__PseudoCaloMETCut", labelData));
    TH1F *if___HLT_MET105_IsoTrk50___RecoPFMET__PseudoCaloMETCut_DATA = (TH1F*)ifileDATA->Get(Form("%s_if___HLT_MET105_IsoTrk50___RecoPFMET__PseudoCaloMETCut", labelData));
    TH1F *if___orMETtrg___RecoPFMET__PseudoCaloMETCut_DATA = (TH1F*)ifileDATA->Get(Form("%s_if___orMETtrg___RecoPFMET__PseudoCaloMETCut", labelData));


    // setup 
    TH1F *eff_HLT_PFMET120_PseudoCaloMET_MC = (TH1F*)if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET_MC->Clone("eff_HLT_PFMET120_PseudoCaloMET_MC");
    TH1F *eff_HLT_PFHT500_PseudoCaloMET_MC = (TH1F*)if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET_MC->Clone("eff_HLT_PFHT500_PseudoCaloMET_MC");
    TH1F *eff_HLT_PFMETNoMu120_PseudoCaloMET_MC = (TH1F*)if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET_MC->Clone("eff_HLT_PFMETNoMu120_PseudoCaloMET_MC");
    TH1F *eff_HLT_MET105_PseudoCaloMET_MC = (TH1F*)if___HLT_MET105_IsoTrk50___PseudoCaloMET_MC->Clone("eff_HLT_MET105_PseudoCaloMET_MC");
    TH1F *eff_orMETtrg_PseudoCaloMET_MC = (TH1F*)if___orMETtrg___PseudoCaloMET_MC->Clone("eff_orMETtrg_PseudoCaloMET_MC");
    
    TH1F *eff_HLT_PFMET120_RecoPFMET_MC = (TH1F*)if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET_MC->Clone("eff_HLT_PFMET120_RecoPFMET_MC");
    TH1F *eff_HLT_PFHT500_RecoPFMET_MC = (TH1F*)if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET_MC->Clone("eff_HLT_PFHT500_RecoPFMET_MC");
    TH1F *eff_HLT_PFMETNoMu120_RecoPFMET_MC = (TH1F*)if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET_MC->Clone("eff_HLT_PFMETNoMu120_RecoPFMET_MC");
    TH1F *eff_HLT_MET105_RecoPFMET_MC = (TH1F*)if___HLT_MET105_IsoTrk50___RecoPFMET_MC->Clone("eff_HLT_MET105_RecoPFMET_MC");
    TH1F *eff_orMETtrg_RecoPFMET_MC = (TH1F*)if___orMETtrg___RecoPFMET_MC->Clone("eff_orMETtrg_RecoPFMET_MC");
    
    TH1F *eff_HLT_PFMET120_PseudoCaloMET__RecoPFMETCut_MC = (TH1F*)if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET__RecoPFMETCut_MC->Clone("eff_HLT_PFMET120_PseudoCaloMET__RecoPFMETCut_MC");
    TH1F *eff_HLT_PFHT500_PseudoCaloMET__RecoPFMETCut_MC = (TH1F*)if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET__RecoPFMETCut_MC->Clone("eff_HLT_PFHT500_PseudoCaloMET__RecoPFMETCut_MC");
    TH1F *eff_HLT_PFMETNoMu120_PseudoCaloMET__RecoPFMETCut_MC = (TH1F*)if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET__RecoPFMETCut_MC->Clone("eff_HLT_PFMETNoMu120_PseudoCaloMET__RecoPFMETCut_MC");
    TH1F *eff_HLT_MET105_PseudoCaloMET__RecoPFMETCut_MC = (TH1F*)if___HLT_MET105_IsoTrk50___PseudoCaloMET__RecoPFMETCut_MC->Clone("eff_HLT_MET105_PseudoCaloMET__RecoPFMETCut_MC");
    TH1F *eff_orMETtrg_PseudoCaloMET__RecoPFMETCut_MC = (TH1F*)if___orMETtrg___PseudoCaloMET__RecoPFMETCut_MC->Clone("eff_orMETtrg_PseudoCaloMET__RecoPFMETCut_MC");

    TH1F *eff_HLT_PFMET120_RecoPFMET__PseudoCaloMETCut_MC = (TH1F*)if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET__PseudoCaloMETCut_MC->Clone("eff_HLT_PFMET120_RecoPFMET__PseudoCaloMETCut_MC");
    TH1F *eff_HLT_PFHT500_RecoPFMET__PseudoCaloMETCut_MC = (TH1F*)if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET__PseudoCaloMETCut_MC->Clone("eff_HLT_PFHT500_RecoPFMET__PseudoCaloMETCut_MC");
    TH1F *eff_HLT_PFMETNoMu120_RecoPFMET__PseudoCaloMETCut_MC = (TH1F*)if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET__PseudoCaloMETCut_MC->Clone("eff_HLT_PFMETNoMu120_RecoPFMET__PseudoCaloMETCut_MC");
    TH1F *eff_HLT_MET105_RecoPFMET__PseudoCaloMETCut_MC = (TH1F*)if___HLT_MET105_IsoTrk50___RecoPFMET__PseudoCaloMETCut_MC->Clone("eff_HLT_MET105_RecoPFMET__PseudoCaloMETCut_MC");
    TH1F *eff_orMETtrg_RecoPFMET__PseudoCaloMETCut_MC = (TH1F*)if___orMETtrg___RecoPFMET__PseudoCaloMETCut_MC->Clone("eff_orMETtrg_RecoPFMET__PseudoCaloMETCut_MC");


    TH1F *eff_HLT_PFMET120_PseudoCaloMET_DATA = (TH1F*)if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET_DATA->Clone("eff_HLT_PFMET120_PseudoCaloMET_DATA");
    TH1F *eff_HLT_PFHT500_PseudoCaloMET_DATA = (TH1F*)if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET_DATA->Clone("eff_HLT_PFHT500_PseudoCaloMET_DATA");
    TH1F *eff_HLT_PFMETNoMu120_PseudoCaloMET_DATA = (TH1F*)if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET_DATA->Clone("eff_HLT_PFMETNoMu120_PseudoCaloMET_DATA");
    TH1F *eff_HLT_MET105_PseudoCaloMET_DATA = (TH1F*)if___HLT_MET105_IsoTrk50___PseudoCaloMET_DATA->Clone("eff_HLT_MET105_PseudoCaloMET_DATA");
    TH1F *eff_orMETtrg_PseudoCaloMET_DATA = (TH1F*)if___orMETtrg___PseudoCaloMET_DATA->Clone("eff_orMETtrg_PseudoCaloMET_DATA");
    
    TH1F *eff_HLT_PFMET120_RecoPFMET_DATA = (TH1F*)if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET_DATA->Clone("eff_HLT_PFMET120_RecoPFMET_DATA");
    TH1F *eff_HLT_PFHT500_RecoPFMET_DATA = (TH1F*)if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET_DATA->Clone("eff_HLT_PFHT500_RecoPFMET_DATA");
    TH1F *eff_HLT_PFMETNoMu120_RecoPFMET_DATA = (TH1F*)if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET_DATA->Clone("eff_HLT_PFMETNoMu120_RecoPFMET_DATA");
    TH1F *eff_HLT_MET105_RecoPFMET_DATA = (TH1F*)if___HLT_MET105_IsoTrk50___RecoPFMET_DATA->Clone("eff_HLT_MET105_RecoPFMET_DATA");
    TH1F *eff_orMETtrg_RecoPFMET_DATA = (TH1F*)if___orMETtrg___RecoPFMET_DATA->Clone("eff_orMETtrg_RecoPFMET_DATA");
    
    TH1F *eff_HLT_PFMET120_PseudoCaloMET__RecoPFMETCut_DATA = (TH1F*)if___HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET__RecoPFMETCut_DATA->Clone("eff_HLT_PFMET120_PseudoCaloMET__RecoPFMETCut_DATA");
    TH1F *eff_HLT_PFHT500_PseudoCaloMET__RecoPFMETCut_DATA = (TH1F*)if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET__RecoPFMETCut_DATA->Clone("eff_HLT_PFHT500_PseudoCaloMET__RecoPFMETCut_DATA");
    TH1F *eff_HLT_PFMETNoMu120_PseudoCaloMET__RecoPFMETCut_DATA = (TH1F*)if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET__RecoPFMETCut_DATA->Clone("eff_HLT_PFMETNoMu120_PseudoCaloMET__RecoPFMETCut_DATA");
    TH1F *eff_HLT_MET105_PseudoCaloMET__RecoPFMETCut_DATA = (TH1F*)if___HLT_MET105_IsoTrk50___PseudoCaloMET__RecoPFMETCut_DATA->Clone("eff_HLT_MET105_PseudoCaloMET__RecoPFMETCut_DATA");
    TH1F *eff_orMETtrg_PseudoCaloMET__RecoPFMETCut_DATA = (TH1F*)if___orMETtrg___PseudoCaloMET__RecoPFMETCut_DATA->Clone("eff_orMETtrg_PseudoCaloMET__RecoPFMETCut_DATA");

    TH1F *eff_HLT_PFMET120_RecoPFMET__PseudoCaloMETCut_DATA = (TH1F*)if___HLT_PFMET120_PFMHT120_IDTight___RecoPFMET__PseudoCaloMETCut_DATA->Clone("eff_HLT_PFMET120_RecoPFMET__PseudoCaloMETCut_DATA");
    TH1F *eff_HLT_PFHT500_RecoPFMET__PseudoCaloMETCut_DATA = (TH1F*)if___HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET__PseudoCaloMETCut_DATA->Clone("eff_HLT_PFHT500_RecoPFMET__PseudoCaloMETCut_DATA");
    TH1F *eff_HLT_PFMETNoMu120_RecoPFMET__PseudoCaloMETCut_DATA = (TH1F*)if___HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET__PseudoCaloMETCut_DATA->Clone("eff_HLT_PFMETNoMu120_RecoPFMET__PseudoCaloMETCut_DATA");
    TH1F *eff_HLT_MET105_RecoPFMET__PseudoCaloMETCut_DATA = (TH1F*)if___HLT_MET105_IsoTrk50___RecoPFMET__PseudoCaloMETCut_DATA->Clone("eff_HLT_MET105_RecoPFMET__PseudoCaloMETCut_DATA");
    TH1F *eff_orMETtrg_RecoPFMET__PseudoCaloMETCut_DATA = (TH1F*)if___orMETtrg___RecoPFMET__PseudoCaloMETCut_DATA->Clone("eff_orMETtrg_RecoPFMET__PseudoCaloMETCut_DATA");



    // efficiency calculation
    eff_HLT_PFMET120_PseudoCaloMET_MC->Divide(PseudoCaloMET_MC);
    eff_HLT_PFHT500_PseudoCaloMET_MC->Divide(PseudoCaloMET_MC);
    eff_HLT_PFMETNoMu120_PseudoCaloMET_MC->Divide(PseudoCaloMET_MC);
    eff_HLT_MET105_PseudoCaloMET_MC->Divide(PseudoCaloMET_MC);
    eff_orMETtrg_PseudoCaloMET_MC->Divide(PseudoCaloMET_MC);
    
    eff_HLT_PFMET120_RecoPFMET_MC->Divide(RecoPFMET_MC);
    eff_HLT_PFHT500_RecoPFMET_MC->Divide(RecoPFMET_MC);
    eff_HLT_PFMETNoMu120_RecoPFMET_MC->Divide(RecoPFMET_MC);
    eff_HLT_MET105_RecoPFMET_MC->Divide(RecoPFMET_MC);
    eff_orMETtrg_RecoPFMET_MC->Divide(RecoPFMET_MC);

    eff_HLT_PFMET120_PseudoCaloMET__RecoPFMETCut_MC->Divide(PseudoCaloMET__RecoPFMETCut_MC);
    eff_HLT_PFHT500_PseudoCaloMET__RecoPFMETCut_MC->Divide(PseudoCaloMET__RecoPFMETCut_MC);
    eff_HLT_PFMETNoMu120_PseudoCaloMET__RecoPFMETCut_MC->Divide(PseudoCaloMET__RecoPFMETCut_MC);
    eff_HLT_MET105_PseudoCaloMET__RecoPFMETCut_MC->Divide(PseudoCaloMET__RecoPFMETCut_MC);
    eff_orMETtrg_PseudoCaloMET__RecoPFMETCut_MC->Divide(PseudoCaloMET__RecoPFMETCut_MC);
    
    eff_HLT_PFMET120_RecoPFMET__PseudoCaloMETCut_MC->Divide(RecoPFMET__PseudoCaloMETCut_MC);
    eff_HLT_PFHT500_RecoPFMET__PseudoCaloMETCut_MC->Divide(RecoPFMET__PseudoCaloMETCut_MC);
    eff_HLT_PFMETNoMu120_RecoPFMET__PseudoCaloMETCut_MC->Divide(RecoPFMET__PseudoCaloMETCut_MC);
    eff_HLT_MET105_RecoPFMET__PseudoCaloMETCut_MC->Divide(RecoPFMET__PseudoCaloMETCut_MC);
    eff_orMETtrg_RecoPFMET__PseudoCaloMETCut_MC->Divide(RecoPFMET__PseudoCaloMETCut_MC);


    eff_HLT_PFMET120_PseudoCaloMET_DATA->Divide(PseudoCaloMET_DATA);
    eff_HLT_PFHT500_PseudoCaloMET_DATA->Divide(PseudoCaloMET_DATA);
    eff_HLT_PFMETNoMu120_PseudoCaloMET_DATA->Divide(PseudoCaloMET_DATA);
    eff_HLT_MET105_PseudoCaloMET_DATA->Divide(PseudoCaloMET_DATA);
    eff_orMETtrg_PseudoCaloMET_DATA->Divide(PseudoCaloMET_DATA);
    
    eff_HLT_PFMET120_RecoPFMET_DATA->Divide(RecoPFMET_DATA);
    eff_HLT_PFHT500_RecoPFMET_DATA->Divide(RecoPFMET_DATA);
    eff_HLT_PFMETNoMu120_RecoPFMET_DATA->Divide(RecoPFMET_DATA);
    eff_HLT_MET105_RecoPFMET_DATA->Divide(RecoPFMET_DATA);
    eff_orMETtrg_RecoPFMET_DATA->Divide(RecoPFMET_DATA);

    eff_HLT_PFMET120_PseudoCaloMET__RecoPFMETCut_DATA->Divide(PseudoCaloMET__RecoPFMETCut_DATA);
    eff_HLT_PFHT500_PseudoCaloMET__RecoPFMETCut_DATA->Divide(PseudoCaloMET__RecoPFMETCut_DATA);
    eff_HLT_PFMETNoMu120_PseudoCaloMET__RecoPFMETCut_DATA->Divide(PseudoCaloMET__RecoPFMETCut_DATA);
    eff_HLT_MET105_PseudoCaloMET__RecoPFMETCut_DATA->Divide(PseudoCaloMET__RecoPFMETCut_DATA);
    eff_orMETtrg_PseudoCaloMET__RecoPFMETCut_DATA->Divide(PseudoCaloMET__RecoPFMETCut_DATA);
    
    eff_HLT_PFMET120_RecoPFMET__PseudoCaloMETCut_DATA->Divide(RecoPFMET__PseudoCaloMETCut_DATA);
    eff_HLT_PFHT500_RecoPFMET__PseudoCaloMETCut_DATA->Divide(RecoPFMET__PseudoCaloMETCut_DATA);
    eff_HLT_PFMETNoMu120_RecoPFMET__PseudoCaloMETCut_DATA->Divide(RecoPFMET__PseudoCaloMETCut_DATA);
    eff_HLT_MET105_RecoPFMET__PseudoCaloMETCut_DATA->Divide(RecoPFMET__PseudoCaloMETCut_DATA);
    eff_orMETtrg_RecoPFMET__PseudoCaloMETCut_DATA->Divide(RecoPFMET__PseudoCaloMETCut_DATA);


    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);

    // drawing
    TCanvas *c_HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET = DrawCanvas(eff_HLT_PFMET120_PseudoCaloMET_DATA, eff_HLT_PFMET120_PseudoCaloMET_MC, "c_HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET", "Pseudo MET [GeV]", "eff. HLT_PFMET120_PFMHT120_IDTight",
    "E1", "E1 same", true, "DATA", "MC", "lep", "lep", 0.045, 0, 1200, 0, 1, false);
    c_HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET->cd();
    latex1->Draw();
    TCanvas *c_HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET = DrawCanvas(eff_HLT_PFHT500_PseudoCaloMET_DATA, eff_HLT_PFHT500_PseudoCaloMET_MC, "c_HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET", "Pseudo MET [GeV]", "eff. HLT_PFHT500_PFMET100_PFMHT100_IDTight",
    "E1", "E1 same", true, "DATA", "MC", "lep", "lep", 0.04, 0, 1200, 0, 1, false);
    c_HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET->cd();
    latex1->Draw();
    TCanvas *c_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET = DrawCanvas(eff_HLT_PFMETNoMu120_PseudoCaloMET_DATA, eff_HLT_PFMETNoMu120_PseudoCaloMET_MC, "c_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET", "Pseudo MET [GeV]", "eff. HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",
    "E1", "E1 same", true, "DATA", "MC", "lep", "lep", 0.032, 0, 1200, 0, 1, false);
    c_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET->cd();
    latex1->Draw();
    TCanvas *c_HLT_MET105_IsoTrk50___PseudoCaloMET = DrawCanvas(eff_HLT_MET105_PseudoCaloMET_DATA, eff_HLT_MET105_PseudoCaloMET_MC, "c_HLT_MET105_IsoTrk50___PseudoCaloMET", "Pseudo MET [GeV]", "eff. HLT_MET105_IsoTrk50",
    "E1", "E1 same", true, "DATA", "MC", "lep", "lep", 0.05, 0, 1200, 0, 1, false);
    c_HLT_MET105_IsoTrk50___PseudoCaloMET->cd();
    latex1->Draw();
    TCanvas *c_orMETtrg___PseudoCaloMET = DrawCanvas(eff_orMETtrg_PseudoCaloMET_DATA, eff_orMETtrg_PseudoCaloMET_MC, "c_orMETtrg___PseudoCaloMET", "Pseudo MET [GeV]", "eff orMETtrg",
    "E1", "E1 same", true, "DATA", "MC", "lep", "lep", 0.06, 0, 1200, 0, 1, false);
    c_orMETtrg___PseudoCaloMET->cd();
    latex1->Draw();

    /*TCanvas *c_HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET = DrawCanvas(eff_HLT_PFHT500_RecoPFMET_DATA, eff_HLT_PFHT500_RecoPFMET_MC, "c_HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET", "RecoPFMET [GeV]", "eff. HLT_PFHT500_PFMET100_PFMHT100_IDTight",
    "E1", "E1 same", true, "DATA", "MC", "lep", "lep", 0, 1200, 0, 1, false);
    TCanvas *c_HLT_PFMET120_PFMHT120_IDTight___RecoPFMET = DrawCanvas(eff_HLT_PFMET120_RecoPFMET_DATA, eff_HLT_PFMET120_RecoPFMET_MC, "c_HLT_PFMET120_PFMHT120_IDTight___RecoPFMET", "RecoPFMET [GeV]", "eff. HLT_PFMET120_PFMHT120_IDTight",
    "E1", "E1 same", true, "DATA", "MC", "lep", "lep", 0, 1200, 0, 1, false);
    TCanvas *c_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET = DrawCanvas(eff_HLT_PFMETNoMu120_RecoPFMET_DATA, eff_HLT_PFMETNoMu120_RecoPFMET_MC, "c_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET", "RecoPFMET [GeV]", "eff. HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",
    "E1", "E1 same", true, "DATA", "MC", "lep", "lep", 0, 1200, 0, 1, false);
    TCanvas *c_HLT_MET105_IsoTrk50___RecoPFMET = DrawCanvas(eff_HLT_MET105_RecoPFMET_DATA, eff_HLT_MET105_RecoPFMET_MC, "c_HLT_MET105_IsoTrk50___RecoPFMET", "RecoPFMET [GeV]", "eff. HLT_MET105_IsoTrk50",
    "E1", "E1 same", true, "DATA", "MC", "lep", "lep", 0, 1200, 0, 1, false);
    TCanvas *c_orMETtrg___RecoPFMET = DrawCanvas(eff_orMETtrg_RecoPFMET_DATA, eff_orMETtrg_RecoPFMET_MC, "c_orMETtrg___RecoPFMET", "RecoPFMET [GeV]", "eff orMETtrg",
    "E1", "E1 same", true, "DATA", "MC", "lep", "lep", 0, 1200, 0, 1, false);

    TCanvas *c_HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET__RecoPFMETCut = DrawCanvas(eff_HLT_PFMET120_PseudoCaloMET__RecoPFMETCut_DATA, eff_HLT_PFMET120_PseudoCaloMET__RecoPFMETCut_MC, "c_HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET__RecoPFMETCut", "Pseudo MET [GeV]", "eff. HLT_PFMET120_PFMHT120_IDTight w/RecoPFMET>170GeV",
    "E1", "E1 same", true, "DATA", "MC", "lep", "lep", 0, 1200, 0, 1, false);
    TCanvas *c_HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET__RecoPFMETCut = DrawCanvas(eff_HLT_PFHT500_PseudoCaloMET__RecoPFMETCut_DATA, eff_HLT_PFHT500_PseudoCaloMET__RecoPFMETCut_MC, "c_HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET__RecoPFMETCut", "Pseudo MET [GeV]", "eff. HLT_PFHT500_PFMET100_PFMHT100_IDTight w/RecoPFMET>170GeV",
    "E1", "E1 same", true, "DATA", "MC", "lep", "lep", 0, 1200, 0, 1, false);
    TCanvas *c_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET__RecoPFMETCut = DrawCanvas(eff_HLT_PFMETNoMu120_PseudoCaloMET__RecoPFMETCut_DATA, eff_HLT_PFMETNoMu120_PseudoCaloMET__RecoPFMETCut_MC, "c_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET__RecoPFMETCut", "Pseudo MET [GeV]", "eff. HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 w/RecoPFMET>170GeV",
    "E1", "E1 same", true, "DATA", "MC", "lep", "lep", 0, 1200, 0, 1, false);
    TCanvas *c_HLT_MET105_IsoTrk50___PseudoCaloMET__RecoPFMETCut = DrawCanvas(eff_HLT_MET105_PseudoCaloMET__RecoPFMETCut_DATA, eff_HLT_MET105_PseudoCaloMET__RecoPFMETCut_MC, "c_HLT_MET105_IsoTrk50___PseudoCaloMET__RecoPFMETCut", "Pseudo MET [GeV]", "eff. HLT_MET105_IsoTrk50 w/RecoPFMET>170GeV",
    "E1", "E1 same", true, "DATA", "MC", "lep", "lep", 0, 1200, 0, 1, false);
    TCanvas *c_orMETtrg___PseudoCaloMET__RecoPFMETCut = DrawCanvas(eff_orMETtrg_PseudoCaloMET__RecoPFMETCut_DATA, eff_orMETtrg_PseudoCaloMET__RecoPFMETCut_MC, "c_orMETtrg___PseudoCaloMET__RecoPFMETCut", "Pseudo MET [GeV]", "eff orMETtrg w/RecoPFMET>170GeV",
    "E1", "E1 same", true, "DATA", "MC", "lep", "lep", 0, 1200, 0, 1, false);

    TCanvas *c_HLT_PFMET120_PFMHT120_IDTight___RecoPFMET__PseudoCaloMETCut = DrawCanvas(eff_HLT_PFMET120_RecoPFMET__PseudoCaloMETCut_DATA, eff_HLT_PFMET120_RecoPFMET__PseudoCaloMETCut_MC, "c_HLT_PFMET120_PFMHT120_IDTight___RecoPFMET__PseudoCaloMETCut", "RecoPFMET [GeV]", "eff. HLT_PFMET120_PFMHT120_IDTight w/PseudoCaloMET>170GeV",
    "E1", "E1 same", true, "DATA", "MC", "lep", "lep", 0, 1200, 0, 1, false);
    TCanvas *c_HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET__PseudoCaloMETCut = DrawCanvas(eff_HLT_PFHT500_RecoPFMET__PseudoCaloMETCut_DATA, eff_HLT_PFHT500_RecoPFMET__PseudoCaloMETCut_MC, "c_HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET__PseudoCaloMETCut", "RecoPFMET [GeV]", "eff. HLT_PFHT500_PFMET100_PFMHT100_IDTight w/PseudoCaloMET>170GeV",
    "E1", "E1 same", true, "DATA", "MC", "lep", "lep", 0, 1200, 0, 1, false);
    TCanvas *c_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET__PseudoCaloMETCut = DrawCanvas(eff_HLT_PFMETNoMu120_RecoPFMET__PseudoCaloMETCut_DATA, eff_HLT_PFMETNoMu120_RecoPFMET__PseudoCaloMETCut_MC, "c_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET__PseudoCaloMETCut", "RecoPFMET [GeV]", "eff. HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 w/PseudoCaloMET>170GeV",
    "E1", "E1 same", true, "DATA", "MC", "lep", "lep", 0, 1200, 0, 1, false);
    TCanvas *c_HLT_MET105_IsoTrk50___RecoPFMET__PseudoCaloMETCut = DrawCanvas(eff_HLT_MET105_RecoPFMET__PseudoCaloMETCut_DATA, eff_HLT_MET105_RecoPFMET__PseudoCaloMETCut_MC, "c_HLT_MET105_IsoTrk50___RecoPFMET__PseudoCaloMETCut", "RecoPFMET [GeV]", "eff. HLT_MET105_IsoTrk50 w/PseudoCaloMET>170GeV",
    "E1", "E1 same", true, "DATA", "MC", "lep", "lep", 0, 1200, 0, 1, false);
    TCanvas *c_orMETtrg___RecoPFMET__PseudoCaloMETCut = DrawCanvas(eff_orMETtrg_RecoPFMET__PseudoCaloMETCut_DATA, eff_orMETtrg_RecoPFMET__PseudoCaloMETCut_MC, "c_orMETtrg___RecoPFMET__PseudoCaloMETCut", "RecoPFMET [GeV]", "eff. orMETtrg w/PseudoCaloMET>170GeV",
    "E1", "E1 same", true, "DATA", "MC", "lep", "lep", 0, 1200, 0, 1, false);*/

    TCanvas *cRatio_HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET = DrawWithRatio(eff_HLT_PFMET120_PseudoCaloMET_DATA, eff_HLT_PFMET120_PseudoCaloMET_MC, c_HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET, "HLT_PFMET120_PFMHT120_IDTight_pseudoCaloMET", "DATA/MC", "Pseudo MET [GeV]", 0, 1200);
    TCanvas *cRatio_HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET = DrawWithRatio(eff_HLT_PFHT500_PseudoCaloMET_DATA, eff_HLT_PFHT500_PseudoCaloMET_MC, c_HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET, "HLT_PFHT500_PFMET100_PFMHT100_IDTight_pseudoCaloMET", "DATA/MC", "Pseudo MET [GeV]", 0, 1200);
    TCanvas *cRatio_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET = DrawWithRatio(eff_HLT_PFMETNoMu120_PseudoCaloMET_DATA, eff_HLT_PFMETNoMu120_PseudoCaloMET_MC, c_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET, "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_pseudoCaloMET", "DATA/MC", "Pseudo MET [GeV]", 0, 1200);
    TCanvas *cRatio_HLT_MET105_IsoTrk50___PseudoCaloMET = DrawWithRatio(eff_HLT_MET105_PseudoCaloMET_DATA, eff_HLT_MET105_PseudoCaloMET_MC, c_HLT_MET105_IsoTrk50___PseudoCaloMET, "HLT_MET105_IsoTrk50_pseudoCaloMET", "DATA/MC", "Pseudo MET [GeV]", 0, 1200);
    TCanvas *cRatio_orMETtrg___PseudoCaloMET = DrawWithRatio(eff_orMETtrg_PseudoCaloMET_DATA, eff_orMETtrg_PseudoCaloMET_MC, c_orMETtrg___PseudoCaloMET, "orMETtrg_pseudoCaloMET", "DATA/MC", "Pseudo MET [GeV]", 0, 1200);
    
    latex1->SetTitle("#it{Private work (CMS simulation/data)}");

    cRatio_HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET->SaveAs(Form("../outputDisplay/cRatio_HLT_PFMET120_PFMHT120_IDTight_%s.pdf", ofilename));
    cRatio_HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET->SaveAs(Form("../outputDisplay/cRatio_HLT_PFHT500_PFMET100_PFMHT100_IDTight_%s.pdf", ofilename));
    cRatio_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET->SaveAs(Form("../outputDisplay/cRatio_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_%s.pdf", ofilename));
    cRatio_HLT_MET105_IsoTrk50___PseudoCaloMET->SaveAs(Form("../outputDisplay/cRatio_HLT_MET105_IsoTrk50_%s.pdf", ofilename));
    cRatio_orMETtrg___PseudoCaloMET->SaveAs(Form("../outputDisplay/cRatio_orMETtrg_%s.pdf", ofilename));


    cRatio_HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET->Modified();
    cRatio_HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET->Update();
    TCanvas *cRatio_HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET_bis = DrawWithRatio(eff_HLT_PFMET120_PseudoCaloMET_DATA, eff_HLT_PFMET120_PseudoCaloMET_MC, c_HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET, "HLT_PFMET120_PFMHT120_IDTight_pseudoCaloMET", "DATA/MC", "Pseudo MET [GeV]", 0, 1200);
    cRatio_HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET_bis->SaveAs(Form("../outputDisplay/cRatio_HLT_PFMET120_PFMHT120_IDTight_%s_bis.pdf", ofilename));
    cRatio_HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET->Modified();
    cRatio_HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET->Update();
    TCanvas *cRatio_HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET_bis = DrawWithRatio(eff_HLT_PFHT500_PseudoCaloMET_DATA, eff_HLT_PFHT500_PseudoCaloMET_MC, c_HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET, "HLT_PFHT500_PFMET100_PFMHT100_IDTight_pseudoCaloMET", "DATA/MC", "Pseudo MET [GeV]", 0, 1200);
    cRatio_HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET_bis->SaveAs(Form("../outputDisplay/cRatio_HLT_PFHT500_PFMET100_PFMHT100_IDTight_%s_bis.pdf", ofilename));
    cRatio_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET->Modified();
    cRatio_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET->Update();
    TCanvas *cRatio_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET_bis = DrawWithRatio(eff_HLT_PFMETNoMu120_PseudoCaloMET_DATA, eff_HLT_PFMETNoMu120_PseudoCaloMET_MC, c_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET, "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_pseudoCaloMET", "DATA/MC", "Pseudo MET [GeV]", 0, 1200);
    cRatio_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET_bis->SaveAs(Form("../outputDisplay/cRatio_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_%s_bis.pdf", ofilename));
    cRatio_HLT_MET105_IsoTrk50___PseudoCaloMET->Modified();
    cRatio_HLT_MET105_IsoTrk50___PseudoCaloMET->Update();
    TCanvas *cRatio_HLT_MET105_IsoTrk50___PseudoCaloMET_bis = DrawWithRatio(eff_HLT_MET105_PseudoCaloMET_DATA, eff_HLT_MET105_PseudoCaloMET_MC, c_HLT_MET105_IsoTrk50___PseudoCaloMET, "HLT_MET105_IsoTrk50_pseudoCaloMET", "DATA/MC", "Pseudo MET [GeV]", 0, 1200);
    cRatio_HLT_MET105_IsoTrk50___PseudoCaloMET_bis->SaveAs(Form("../outputDisplay/cRatio_HLT_MET105_IsoTrk50_%s_bis.pdf", ofilename));
    cRatio_orMETtrg___PseudoCaloMET->Modified();
    cRatio_orMETtrg___PseudoCaloMET->Update();
    TCanvas *cRatio_orMETtrg___PseudoCaloMET_bis = DrawWithRatio(eff_orMETtrg_PseudoCaloMET_DATA, eff_orMETtrg_PseudoCaloMET_MC, c_orMETtrg___PseudoCaloMET, "orMETtrg_pseudoCaloMET", "DATA/MC", "Pseudo MET [GeV]", 0, 1200);
    cRatio_orMETtrg___PseudoCaloMET_bis->SaveAs(Form("../outputDisplay/cRatio_orMETtrg_%s_bis.pdf", ofilename));



    /*TCanvas *cRatio_HLT_PFMET120_PFMHT120_IDTight___RecoPFMET = DrawWithRatio(eff_HLT_PFMET120_RecoPFMET_DATA, eff_HLT_PFMET120_RecoPFMET_MC, c_HLT_PFMET120_PFMHT120_IDTight___RecoPFMET, "HLT_PFMET120_PFMHT120_IDTight_pfMET", "DATA/MC", "RecoPFMET [GeV]", 0, 1200);
    TCanvas *cRatio_HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET = DrawWithRatio(eff_HLT_PFHT500_RecoPFMET_DATA, eff_HLT_PFHT500_RecoPFMET_MC, c_HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET, "HLT_PFHT500_PFMET100_PFMHT100_IDTight_pfMET", "DATA/MC", "RecoPFMET [GeV]", 0, 1200);
    TCanvas *cRatio_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET = DrawWithRatio(eff_HLT_PFMETNoMu120_RecoPFMET_DATA, eff_HLT_PFMETNoMu120_RecoPFMET_MC, c_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET, "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_pfMET", "DATA/MC", "RecoPFMET [GeV]", 0, 1200);
    TCanvas *cRatio_HLT_MET105_IsoTrk50___RecoPFMET = DrawWithRatio(eff_HLT_MET105_RecoPFMET_DATA, eff_HLT_MET105_RecoPFMET_MC, c_HLT_MET105_IsoTrk50___RecoPFMET, "HLT_MET105_IsoTrk_pfMET", "DATA/MC", "RecoPFMET [GeV]", 0, 1200);
    TCanvas *cRatio_orMETtrg___RecoPFMET = DrawWithRatio(eff_orMETtrg_RecoPFMET_DATA, eff_orMETtrg_RecoPFMET_MC, c_orMETtrg___RecoPFMET, "orMETtrg_pfMET", "DATA/MC", "RecoPFMET [GeV]", 0, 1200);
    
    TCanvas *cRatio_HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET__RecoPFMETCut = DrawWithRatio(eff_HLT_PFMET120_PseudoCaloMET__RecoPFMETCut_DATA, eff_HLT_PFMET120_PseudoCaloMET__RecoPFMETCut_MC, c_HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET__RecoPFMETCut, "HLT_PFMET120_PFMHT120_IDTight_w/RecoPFMET>170GeV", "DATA/MC", "Pseudo MET [GeV]", 0, 1200);
    TCanvas *cRatio_HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET__RecoPFMETCut = DrawWithRatio(eff_HLT_PFHT500_PseudoCaloMET__RecoPFMETCut_DATA, eff_HLT_PFHT500_PseudoCaloMET__RecoPFMETCut_MC, c_HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET__RecoPFMETCut, "HLT_PFHT500_PFMET100_PFMHT100_IDTight_w/RecoPFMET>170GeV", "DATA/MC", "Pseudo MET [GeV]", 0, 1200);
    TCanvas *cRatio_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET__RecoPFMETCut = DrawWithRatio(eff_HLT_PFMETNoMu120_PseudoCaloMET__RecoPFMETCut_DATA, eff_HLT_PFMETNoMu120_PseudoCaloMET__RecoPFMETCut_MC, c_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET__RecoPFMETCut, "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_w/RecoPFMET>170GeV", "DATA/MC", "Pseudo MET [GeV]", 0, 1200);
    TCanvas *cRatio_HLT_MET105_IsoTrk50___PseudoCaloMET__RecoPFMETCut = DrawWithRatio(eff_HLT_MET105_PseudoCaloMET__RecoPFMETCut_DATA, eff_HLT_MET105_PseudoCaloMET__RecoPFMETCut_MC, c_HLT_MET105_IsoTrk50___PseudoCaloMET__RecoPFMETCut, "HLT_MET105_IsoTrk50_w/RecoPFMET>170GeV", "DATA/MC", "Pseudo MET [GeV]", 0, 1200);
    TCanvas *cRatio_orMETtrg___PseudoCaloMET__RecoPFMETCut = DrawWithRatio(eff_orMETtrg_PseudoCaloMET__RecoPFMETCut_DATA, eff_orMETtrg_PseudoCaloMET__RecoPFMETCut_MC, c_orMETtrg___PseudoCaloMET__RecoPFMETCut, "orMETtrg_w/RecoPFMET>170GeV", "DATA/MC", "Pseudo MET [GeV]", 0, 1200);

    TCanvas *cRatio_HLT_PFMET120_PFMHT120_IDTight___RecoPFMET__PseudoCaloMETCut = DrawWithRatio(eff_HLT_PFMET120_RecoPFMET__PseudoCaloMETCut_DATA, eff_HLT_PFMET120_RecoPFMET__PseudoCaloMETCut_MC, c_HLT_PFMET120_PFMHT120_IDTight___RecoPFMET__PseudoCaloMETCut, "HLT_PFMET120_PFMHT120_IDTight_w/PseudoCaloMET>170GeV", "DATA/MC", "RecoPFMET [GeV]", 0, 1200);
    TCanvas *cRatio_HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET__PseudoCaloMETCut = DrawWithRatio(eff_HLT_PFHT500_RecoPFMET__PseudoCaloMETCut_DATA, eff_HLT_PFHT500_RecoPFMET__PseudoCaloMETCut_MC, c_HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET__PseudoCaloMETCut, "HLT_PFHT500_PFMET100_PFMHT100_IDTight_w/PseudoCaloMET>170GeV", "DATA/MC", "RecoPFMET [GeV]", 0, 1200);
    TCanvas *cRatio_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET__PseudoCaloMETCut = DrawWithRatio(eff_HLT_PFMETNoMu120_RecoPFMET__PseudoCaloMETCut_DATA, eff_HLT_PFMETNoMu120_RecoPFMET__PseudoCaloMETCut_MC, c_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET__PseudoCaloMETCut, "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_w/PseudoCaloMET>170GeV", "DATA/MC", "RecoPFMET [GeV]", 0, 1200);
    TCanvas *cRatio_HLT_MET105_IsoTrk50___RecoPFMET__PseudoCaloMETCut = DrawWithRatio(eff_HLT_MET105_RecoPFMET__PseudoCaloMETCut_DATA, eff_HLT_MET105_RecoPFMET__PseudoCaloMETCut_MC, c_HLT_MET105_IsoTrk50___RecoPFMET__PseudoCaloMETCut, "HLT_MET105_IsoTrk50_w/PseudoCaloMET>170GeV", "DATA/MC", "RecoPFMET [GeV]", 0, 1200);
    TCanvas *cRatio_orMETtrg___RecoPFMET__PseudoCaloMETCut = DrawWithRatio(eff_orMETtrg_RecoPFMET__PseudoCaloMETCut_DATA, eff_orMETtrg_RecoPFMET__PseudoCaloMETCut_MC, c_orMETtrg___RecoPFMET__PseudoCaloMETCut, "orMETtrg_w/PseudoCaloMET>170GeV", "DATA/MC", "RecoPFMET [GeV]", 0, 1200);
    

    // saving
    ofile->cd();
    cRatio_HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET->Write();
    cRatio_HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET->Write();
    cRatio_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET->Write();
    cRatio_HLT_MET105_IsoTrk50___PseudoCaloMET->Write();
    cRatio_orMETtrg___PseudoCaloMET->Write();
    
    cRatio_HLT_PFMET120_PFMHT120_IDTight___RecoPFMET->Write();
    cRatio_HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET->Write();
    cRatio_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET->Write();
    cRatio_HLT_MET105_IsoTrk50___RecoPFMET->Write();
    cRatio_orMETtrg___RecoPFMET->Write();
    
    cRatio_HLT_PFMET120_PFMHT120_IDTight___PseudoCaloMET__RecoPFMETCut->Write();
    cRatio_HLT_PFHT500_PFMET100_PFMHT100_IDTight___PseudoCaloMET__RecoPFMETCut->Write();
    cRatio_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___PseudoCaloMET__RecoPFMETCut->Write();
    cRatio_HLT_MET105_IsoTrk50___PseudoCaloMET__RecoPFMETCut->Write();
    cRatio_orMETtrg___PseudoCaloMET__RecoPFMETCut->Write();

    cRatio_HLT_PFMET120_PFMHT120_IDTight___RecoPFMET__PseudoCaloMETCut->Write();
    cRatio_HLT_PFHT500_PFMET100_PFMHT100_IDTight___RecoPFMET__PseudoCaloMETCut->Write();
    cRatio_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60___RecoPFMET__PseudoCaloMETCut->Write();
    cRatio_HLT_MET105_IsoTrk50___RecoPFMET__PseudoCaloMETCut->Write();
    cRatio_orMETtrg___RecoPFMET__PseudoCaloMETCut->Write();
    ofile->Close();*/


    return;
}


void Comp_muonEG(const char *inputfileDATA, const char *inputfileMC) {

    TFile *ofile = new TFile("../outputDisplay/Comp_muonEG.root", "RECREATE");

    TFile *ifileDATA = new TFile(Form("%s", inputfileDATA), "READ");
    TFile *ifileMC = new TFile(Form("%s", inputfileMC), "READ");
    // histograms
    TH1F *electron_pt_DATA = (TH1F*)ifileDATA->Get("CalibPseudoMET_electron_pt");
    TH1F *electron_eta_DATA = (TH1F*)ifileDATA->Get("CalibPseudoMET_electron_eta");
    TH1F *electron_phi_DATA = (TH1F*)ifileDATA->Get("CalibPseudoMET_electron_phi");
    TH1F *muon_pt_DATA = (TH1F*)ifileDATA->Get("CalibPseudoMET_muon_pt");
    TH1F *muon_eta_DATA = (TH1F*)ifileDATA->Get("CalibPseudoMET_muon_eta");
    TH1F *muon_phi_DATA = (TH1F*)ifileDATA->Get("CalibPseudoMET_muon_phi");    
    TH1F *PseudoCaloMET_DATA = (TH1F*)ifileDATA->Get("CalibPseudoMET_PseudoCaloMET");
    TH1F *RecoPFMET_DATA = (TH1F*)ifileDATA->Get("CalibPseudoMET_RecoPFMET");
    TH1F *L1MET_DATA = (TH1F*)ifileDATA->Get("CalibPseudoMET_L1MET");
    TH1F *HLTCaloMET_DATA = (TH1F*)ifileDATA->Get("CalibPseudoMET_HLTCaloMET");
    TH1F *HLTCaloMHT_DATA = (TH1F*)ifileDATA->Get("CalibPseudoMET_HLTCaloMHT");
    TH1F *HLTPFMHT_DATA = (TH1F*)ifileDATA->Get("CalibPseudoMET_HLTPFMHT");
    TH1F *HLTPFMET_DATA = (TH1F*)ifileDATA->Get("CalibPseudoMET_HLTPFMET");

    TH1F *electron_pt_MC = (TH1F*)ifileMC->Get("CalibPseudoMET_electron_pt");
    TH1F *electron_eta_MC = (TH1F*)ifileMC->Get("CalibPseudoMET_electron_eta");
    TH1F *electron_phi_MC = (TH1F*)ifileMC->Get("CalibPseudoMET_electron_phi");
    TH1F *muon_pt_MC = (TH1F*)ifileMC->Get("CalibPseudoMET_muon_pt");
    TH1F *muon_eta_MC = (TH1F*)ifileMC->Get("CalibPseudoMET_muon_eta");
    TH1F *muon_phi_MC = (TH1F*)ifileMC->Get("CalibPseudoMET_muon_phi");
    TH1F *PseudoCaloMET_MC = (TH1F*)ifileMC->Get("CalibPseudoMET_PseudoCaloMET");
    TH1F *RecoPFMET_MC = (TH1F*)ifileMC->Get("CalibPseudoMET_RecoPFMET");
    TH1F *L1MET_MC = (TH1F*)ifileMC->Get("CalibPseudoMET_L1MET");
    TH1F *HLTCaloMET_MC = (TH1F*)ifileMC->Get("CalibPseudoMET_HLTCaloMET");
    TH1F *HLTCaloMHT_MC = (TH1F*)ifileMC->Get("CalibPseudoMET_HLTCaloMHT");
    TH1F *HLTPFMHT_MC = (TH1F*)ifileMC->Get("CalibPseudoMET_HLTPFMHT");
    TH1F *HLTPFMET_MC = (TH1F*)ifileMC->Get("CalibPseudoMET_HLTPFMET");


    // scale
    electron_pt_DATA->Scale(1./electron_pt_DATA->Integral());
    electron_eta_DATA->Scale(1./electron_eta_DATA->Integral());
    electron_phi_DATA->Scale(1./electron_phi_DATA->Integral());
    muon_pt_DATA->Scale(1./muon_pt_DATA->Integral());
    muon_eta_DATA->Scale(1./muon_eta_DATA->Integral());
    muon_phi_DATA->Scale(1./muon_phi_DATA->Integral());
    PseudoCaloMET_DATA->Scale(1./PseudoCaloMET_DATA->Integral());
    RecoPFMET_DATA->Scale(1./RecoPFMET_DATA->Integral());
    L1MET_DATA->Scale(1./L1MET_DATA->Integral());
    HLTCaloMET_DATA->Scale(1./HLTCaloMET_DATA->Integral());
    HLTCaloMHT_DATA->Scale(1./HLTCaloMHT_DATA->Integral());
    HLTPFMHT_DATA->Scale(1./HLTPFMHT_DATA->Integral());
    HLTPFMET_DATA->Scale(1./HLTPFMET_DATA->Integral());

    electron_pt_MC->Scale(1./electron_pt_MC->Integral());
    electron_eta_MC->Scale(1./electron_eta_MC->Integral());
    electron_phi_MC->Scale(1./electron_phi_MC->Integral());
    muon_pt_MC->Scale(1./muon_pt_MC->Integral());
    muon_eta_MC->Scale(1./muon_eta_MC->Integral());
    muon_phi_MC->Scale(1./muon_phi_MC->Integral());
    PseudoCaloMET_MC->Scale(1./PseudoCaloMET_MC->Integral());
    RecoPFMET_MC->Scale(1./RecoPFMET_MC->Integral());
    L1MET_MC->Scale(1./L1MET_MC->Integral());
    HLTCaloMET_MC->Scale(1./HLTCaloMET_MC->Integral());
    HLTCaloMHT_MC->Scale(1./HLTCaloMHT_MC->Integral());
    HLTPFMHT_MC->Scale(1./HLTPFMHT_MC->Integral());
    HLTPFMET_MC->Scale(1./HLTPFMET_MC->Integral());


    // rebin
    muon_pt_DATA->Rebin(4);
    muon_pt_MC->Rebin(4);
    electron_pt_DATA->Rebin(4);
    electron_pt_MC->Rebin(4);

    // canvas 
    TCanvas *c_electron_pt = DrawCanvas(electron_pt_DATA, electron_pt_MC, "c_electron_pt", "Electron p_{T} [GeV]", "Number of tracks (normalized)",
    "E1", "HIST same", true, "DATA", "MC", "lep", "f", 0, 600, 0, -1, true, true);
    TCanvas *c_electron_eta = DrawCanvas(electron_eta_DATA, electron_eta_MC, "c_electron_eta", "Electron #eta", "Number of tracks (normalized)",
    "E1", "HIST same", true, "DATA", "MC", "lep", "f", -3.2, 3.2, 0, -1, false, true);
    TCanvas *c_electron_phi = DrawCanvas(electron_phi_DATA, electron_phi_MC, "c_electron_phi", "Electron #phi", "Number of tracks (normalized)",
    "E1", "HIST same", true, "DATA", "MC", "lep", "f", -3.2, 3.2, 0, -1, false, true);

    TCanvas *c_muon_pt = DrawCanvas(muon_pt_DATA, muon_pt_MC, "c_muon_pt", "Muon p_{T} [GeV]", "Number of tracks (normalized)",
    "E1", "HIST same", true, "DATA", "MC", "lep", "f", 0, 600, 0, -1, true, true);
    TCanvas *c_muon_eta = DrawCanvas(muon_eta_DATA, muon_eta_MC, "c_muon_eta", "Muon #eta", "Number of tracks (normalized)",
    "E1", "HIST same", true, "DATA", "MC", "lep", "f", -3.2, 3.2, 0, -1, false, true);
    TCanvas *c_muon_phi = DrawCanvas(muon_phi_DATA, muon_phi_MC, "c_muon_phi", "Muon #phi", "Number of tracks (normalized)",
    "E1", "HIST same", true, "DATA", "MC", "lep", "f", -3.2, 3.2, 0, -1, false, true);

    TCanvas *c_PseudoCaloMET = DrawCanvas(PseudoCaloMET_DATA, PseudoCaloMET_MC, "c_PseudoCaloMET", "Pseudo MET [GeV]", "Number of events (normalized)",
    "E1", "HIST same", true, "DATA", "MC", "lep", "f", 0, 1000, 0, -1, true, true);
    TCanvas *c_RecoPFMET = DrawCanvas(RecoPFMET_DATA, RecoPFMET_MC, "c_RecoPFMET", "RecoPFMET [GeV]", "Number of events (normalized)",
    "E1", "HIST same", true, "DATA", "MC", "lep", "f", 0, 1000, 0, -1, true, true);
    TCanvas *c_L1MET = DrawCanvas(L1MET_DATA, L1MET_MC, "c_L1MET", "L1MET [GeV]", "Number of events (normalized)",
    "E1", "HIST same", true, "DATA", "MC", "lep", "f", 0, 1000, 0, -1, true, true);
    TCanvas *c_HLTCaloMET = DrawCanvas(HLTCaloMET_DATA, HLTCaloMET_MC, "c_HLTCaloMET", "HLTCaloMET [GeV]", "Number of events (normalized)",
    "E1", "HIST same", true, "DATA", "MC", "lep", "f", 0, 1000, 0, -1, true, true);
    TCanvas *c_HLTCaloMHT = DrawCanvas(HLTCaloMHT_DATA, HLTCaloMHT_MC, "c_HLTCaloMHT", "HLTCaloMHT [GeV]", "Number of events (normalized)",
    "E1", "HIST same", true, "DATA", "MC", "lep", "f", 0, 1000, 0, -1, true, true);
    TCanvas *c_HLTPFMHT = DrawCanvas(HLTPFMHT_DATA, HLTPFMHT_MC, "c_HLTPFMHT", "HLTPFMHT [GeV]", "Number of events (normalized)",
    "E1", "HIST same", true, "DATA", "MC", "lep", "f", 0, 1000, 0, -1, true, true);
    TCanvas *c_HLTPFMET = DrawCanvas(HLTPFMET_DATA, HLTPFMET_MC, "c_HLTPFMET", "HLTPFMET [GeV]", "Number of events (normalized)",
    "E1", "HIST same", true, "DATA", "MC", "lep", "f", 0, 1000, 0, -1, true, true);


    //draw with ratio
    TCanvas *c_RATIO_electron_pt = DrawWithRatio(electron_pt_DATA, electron_pt_MC, c_electron_pt, "Electron_pT", "DATA/MC", "Electron p_{T}", 0, 600);
    TCanvas *c_RATIO_electron_eta = DrawWithRatio(electron_eta_DATA, electron_eta_MC, c_electron_eta, "Electron_Eta", "DATA/MC", "Electron #eta", -3.2, 3.2);
    TCanvas *c_RATIO_electron_phi = DrawWithRatio(electron_phi_DATA, electron_phi_MC, c_electron_phi, "Electron_Phi", "DATA/MC", "Electron #phi", -3.2, 3.2);

    TCanvas *c_RATIO_muon_pt = DrawWithRatio(muon_pt_DATA, muon_pt_MC, c_muon_pt, "Muon_pT", "DATA/MC", "Muon p_{T}", 0, 600);
    TCanvas *c_RATIO_muon_eta = DrawWithRatio(muon_eta_DATA, muon_eta_MC, c_muon_eta, "Muon_Eta", "DATA/MC", "Muon #eta", -3.2, 3.2);
    TCanvas *c_RATIO_muon_phi = DrawWithRatio(muon_phi_DATA, muon_phi_MC, c_muon_phi, "Muon_Phi", "DATA/MC", "Muon #phi", -3.2, 3.2);

    TCanvas *c_RATIO_PseudoCaloMET = DrawWithRatio(PseudoCaloMET_DATA, PseudoCaloMET_MC, c_PseudoCaloMET, "PseudoCaloMET", "DATA/MC", "Pseudo MET [GeV]", 0, 1000);
    TCanvas *c_RATIO_RecoPFMET = DrawWithRatio(RecoPFMET_DATA, RecoPFMET_MC, c_RecoPFMET, "RecoPFMET", "DATA/MC", "RecoPFMET [GeV]", 0, 1000);
    TCanvas *c_RATIO_L1MET = DrawWithRatio(L1MET_DATA, L1MET_MC, c_L1MET, "L1MET", "DATA/MC", "L1MET [GeV]", 0, 1000);
    TCanvas *c_RATIO_HLTCaloMET = DrawWithRatio(HLTCaloMET_DATA, HLTCaloMET_MC, c_HLTCaloMET, "HLTCaloMET", "DATA/MC", "HLTCaloMET [GeV]", 0, 1000);
    TCanvas *c_RATIO_HLTCaloMHT = DrawWithRatio(HLTCaloMHT_DATA, HLTCaloMHT_MC, c_HLTCaloMHT, "HLTCaloMHT", "DATA/MC", "HLTCaloMHT [GeV]", 0, 1000);
    TCanvas *c_RATIO_HLTPFMHT = DrawWithRatio(HLTPFMHT_DATA, HLTPFMHT_MC, c_HLTPFMHT, "HLTPFMHT", "DATA/MC", "HLTPFMHT [GeV]", 0, 1000);
    TCanvas *c_RATIO_HLTPFMET = DrawWithRatio(HLTPFMET_DATA, HLTPFMET_MC, c_HLTPFMET, "HLTPFMET", "DATA/MC", "HLTPFMET [GeV]", 0, 1000);


    // saving
    ofile->cd();
    c_RATIO_electron_pt->Write();
    c_RATIO_electron_eta->Write();
    c_RATIO_electron_phi->Write();

    c_RATIO_muon_pt->Write();
    c_RATIO_muon_eta->Write();
    c_RATIO_muon_phi->Write();

    c_RATIO_PseudoCaloMET->Write();
    c_RATIO_RecoPFMET->Write();
    c_RATIO_L1MET->Write();
    c_RATIO_HLTCaloMET->Write();
    c_RATIO_HLTCaloMHT->Write();
    c_RATIO_HLTPFMHT->Write();
    c_RATIO_HLTPFMET->Write();
    ofile->Close();

    return;
}

void Old_vs_New_fits(const char *inputfile) {

    TFile *ofile = new TFile("../outputDisplay/Old_vs_New_fits.root", "RECREATE");
    TFile *ifile = new TFile(Form("%s", inputfile), "READ");

    // histograms
    TH1F* Ih = (TH1F*)ifile->Get("METanalysis_Eta2p4_Ih");
    TH1F* oP = (TH1F*)ifile->Get("METanalysis_Eta2p4_10000oP");

    TH1F* Ih_new = (TH1F*)Ih->Clone("Ih_new");
    TH1F* oP_new = (TH1F*)oP->Clone("oP_new");


    // 1oP fit
    float rangemax_p = 30;
    if (oP->GetBinCenter(oP->GetMaximumBin()) < rangemax_p) rangemax_p = 0.8 * oP->GetBinCenter(oP->GetMaximumBin());
    TF1 f_p_old("f_p_old","[0]*([1]+erf((log(x)-[2])/[3]))",0,rangemax_p);
    f_p_old.SetParameter(0,560);
    f_p_old.FixParameter(1,1.0);
    f_p_old.SetParameter(2,3.50116e+00);
    f_p_old.SetParameter(3,0.60152e+00);
    oP->Fit(&f_p_old,"RL","",0,rangemax_p);

    TF1 f_p_new("f_p_new", "0.5*(exp([0]*x*x+[1]*x)+exp(-[0]*x*x-[1]*x))-1", 0, rangemax_p);
    float end1oPFit = 0.4 * oP->GetBinCenter(oP->GetMaximumBin());
    if (end1oPFit > 25) end1oPFit = 25;
    oP_new->Fit(&f_p_new, "R", "", 0, end1oPFit);

    // Ih fit
    float max_ih = Ih_new->GetBinCenter(Ih_new->GetMaximumBin());
    TF1 f_ih_old("f_ih_old", "gaus", 3, 8);
    f_ih_old.SetParameter(0, 0.5*Ih->Integral());
    f_ih_old.SetParameter(1, max_ih);
    f_ih_old.SetParameter(2, Ih->GetStdDev());
    Ih->Fit(&f_ih_old, "RL", "", 3, 8);

    float start_fit = 1.2*max_ih;
    int lastBinContent = Ih->GetNbinsX();
    while(Ih->GetBinContent(lastBinContent)==0) lastBinContent--;
    if(start_fit > Ih->GetBinCenter(lastBinContent)) start_fit = max_ih;
    TF1 f_ih_new("f_ih_new", "gaus", start_fit, 6);
    f_ih_new.SetParameter(0, 0.5*Ih->Integral());
    f_ih_new.SetParameter(1, max_ih);
    f_ih_new.SetParameter(2, Ih->GetStdDev());
    Ih_new->Fit(&f_ih_new, "RL", "", start_fit, 6);


    TCanvas *cIhold = DrawCanvas(Ih, "Ih old fit", "Ih [MeV/cm]", "Number of tracks", "E1", 0, 10, 0, -1, false);
    TCanvas *coPold = DrawCanvas(oP, "oP old fit", "10^{4}/p [GeV^{-1}]", "Number of tracks", "E1", 0, 200, 0, -1, false);
    TCanvas *cIhnew = DrawCanvas(Ih_new, "Ih new fit", "Ih [MeV/cm]", "Number of tracks", "E1", 0, 10, 0, -1, false);
    TCanvas *coPnew = DrawCanvas(oP_new, "oP new fit", "10^{4}/p [GeV^{-1}]", "Number of tracks", "E1", 0, 200, 0, -1, false);

    ofile->cd();
    cIhold->Write();
    coPold->Write();
    cIhnew->Write();
    coPnew->Write();
    ofile->Close();


    return;
}

void Ihand1oP_fits() {

    TFile *ifile = new TFile("../output/JetMET2024_V12/JetMET2024_V12p24.root", "READ");

    // histograms
    TH1F* Ih_new = (TH1F*)ifile->Get("METanalysis_Eta2p4_Ih");
    TH1F* oP_new = (TH1F*)ifile->Get("METanalysis_Eta2p4_10000oP");

    // 1oP fit
    float rangemax_p = 0.9 * oP_new->GetBinCenter(oP_new->GetMaximumBin());
    TF1 f_p_new("f_p_old","[0]*([1]+erf((log(x)-[2])/[3]))",0,rangemax_p);
    f_p_new.SetParameter(0,560);
    f_p_new.FixParameter(1,1.0);
    f_p_new.SetParameter(2,3.50116e+00);
    f_p_new.SetParameter(3,0.60152e+00);
    
    oP_new->Fit(&f_p_new, "RL", "", 0, rangemax_p);

    // Ih fit
    float max_ih = Ih_new->GetBinCenter(Ih_new->GetMaximumBin());
   
    float start_fit = 1.2*max_ih;
    int lastBinContent = Ih_new->GetNbinsX();
    while(Ih_new->GetBinContent(lastBinContent)==0) lastBinContent--;
    if(start_fit > Ih_new->GetBinCenter(lastBinContent)) start_fit = max_ih;
    TF1 f_ih_new("f_ih_new", "gaus", start_fit, 6);
    f_ih_new.SetParameter(0, 0.5*Ih_new->Integral());
    f_ih_new.SetParameter(1, max_ih);
    f_ih_new.SetParameter(2, Ih_new->GetStdDev());
    Ih_new->Fit(&f_ih_new, "RL", "", start_fit, 6);


    Ih_new->GetYaxis()->SetTitleSize(0.06);
    Ih_new->GetXaxis()->SetTitleSize(0.06);
    Ih_new->GetXaxis()->SetTitleOffset(0.9);
    Ih_new->GetYaxis()->SetTitleOffset(1);
    Ih_new->GetXaxis()->SetLabelSize(0.05);
    Ih_new->GetYaxis()->SetLabelSize(0.05);
    Ih_new->GetXaxis()->SetTitle("I_{h} [MeV/cm]");
    Ih_new->GetYaxis()->SetTitle("Number of tracks");
    Ih_new->SetLineColor(kBlack);
    Ih_new->SetMarkerStyle(20);
    Ih_new->SetMarkerColor(kBlack);
    Ih_new->GetXaxis()->SetRangeUser(0, 6);

    oP_new->GetYaxis()->SetTitleSize(0.06);
    oP_new->GetXaxis()->SetTitleSize(0.06);
    oP_new->GetXaxis()->SetTitleOffset(0.9);
    oP_new->GetYaxis()->SetTitleOffset(1);
    oP_new->GetXaxis()->SetLabelSize(0.05);
    oP_new->GetYaxis()->SetLabelSize(0.05);
    oP_new->GetXaxis()->SetTitle("10^{4}/p [GeV^{-1}]");
    oP_new->GetYaxis()->SetTitle("Number of tracks");
    oP_new->SetLineColor(kBlack);
    oP_new->SetMarkerStyle(20);
    oP_new->SetMarkerColor(kBlack);
    oP_new->GetXaxis()->SetRangeUser(0, 210);

    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);

    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    c1->cd();
    c1->SetLeftMargin(0.16); c1->SetBottomMargin(0.16);
    gStyle->SetOptStat(0);
    Ih_new->Draw("EO");
    latex1->Draw();
    c1->SetLogy();
    gStyle->SetOptFit(1);
    f_ih_new.SetLineColor(kRed);
    f_ih_new.SetLineWidth(2);
    f_ih_new.Draw("same");
    
    TPaveStats *st = (TPaveStats*)Ih_new->FindObject("stats");
    if (st) {
        st->SetX1NDC(0.2); // new x start position
        st->SetX2NDC(0.6); // new x end position
        st->SetY1NDC(0.2);  // new y start position
        st->SetY2NDC(0.5);  // new y end position
    }
    st->SetTextSize(0.035);
    c1->SaveAs("../outputDisplay/Ih_newfit.pdf");
    latex1->SetTitle("#it{Private work (CMS data)}");
    c1->Modified();
    c1->Update();
    c1->SaveAs("../outputDisplay/Ih_newfit_bis.pdf");
    latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Work in progress}");


    TCanvas *c2 = new TCanvas("c2","c2",800,600);
    c2->cd();
    c2->SetLeftMargin(0.16); c2->SetBottomMargin(0.16);
    gStyle->SetOptStat(0);
    oP_new->Draw("E0");
    latex1->Draw();
    c2->SetLogy();
    gStyle->SetOptFit(1);
    f_p_new.SetLineColor(kRed);
    f_p_new.SetLineWidth(2);
    f_p_new.SetRange(0, 11);
    f_p_new.Draw("same");
    
    TPaveStats *st2 = (TPaveStats*)oP_new->FindObject("stats");
    if (st2) {
        st2->SetX1NDC(0.3); // new x start position
        st2->SetX2NDC(0.6); // new x end position
        st2->SetY1NDC(0.2);  // new y start position
        st2->SetY2NDC(0.5);  // new y end position
    }
    st2->SetTextSize(0.035);
    c2->SaveAs("../outputDisplay/oP_newfit.pdf");
    latex1->SetTitle("#it{Private work (CMS data)}");
    c2->Modified();
    c2->Update();
    c2->SaveAs("../outputDisplay/oP_newfit_bis.pdf");


    return;
}

void BKGdependency_preliminary (const char *inputname, const char *ofilename) {

    TFile *ofile = new TFile(Form("../outputDisplay/BKGdependency_preliminary_%s.root", ofilename), "RECREATE");
    TFile *ifile = new TFile(Form("%s", inputname), "READ");

    TH2F *pT_vs_fpix = (TH2F*)ifile->Get("METanalysis_Eta2p4_pT_vs_Fpixel");

    TH2F *ih_vs_eta_A3fp9 = (TH2F*)ifile->Get("ih_eta_regionA_3fp9_METanalysis_Eta2p4");
    TH2F *ih_vs_eta_A9fp10 = (TH2F*)ifile->Get("ih_eta_regionA_9fp10_METanalysis_Eta2p4");
    TH2F *ih_vs_eta_D3fp8 = (TH2F*)ifile->Get("ih_eta_regionD_3fp8_METanalysis_Eta2p4");
    TH2F *ih_vs_eta_D8fp9 = (TH2F*)ifile->Get("ih_eta_regionD_8fp9_METanalysis_Eta2p4");
    TH2F *ih_vs_eta_D9fp10 = (TH2F*)ifile->Get("ih_eta_regionD_9fp10_METanalysis_Eta2p4");

    TH2F *fpix_vs_ih_A3fp9 = (TH2F*)ifile->Get("fpix_ih_regionA_3fp9_METanalysis_Eta2p4");
    TH2F *fpix_vs_ih_A9fp10 = (TH2F*)ifile->Get("fpix_ih_regionA_9fp10_METanalysis_Eta2p4");
    TH2F *fpix_vs_ih_D3fp8 = (TH2F*)ifile->Get("fpix_ih_regionD_3fp8_METanalysis_Eta2p4");
    TH2F *fpix_vs_ih_D8fp9 = (TH2F*)ifile->Get("fpix_ih_regionD_8fp9_METanalysis_Eta2p4");
    TH2F *fpix_vs_ih_D9fp10 = (TH2F*)ifile->Get("fpix_ih_regionD_9fp10_METanalysis_Eta2p4");

    TH2F *oP_vs_fpix_A3fp9 = (TH2F*)ifile->Get("oP_fpix_regionA_3fp9_METanalysis_Eta2p4");
    TH2F *oP_vs_fpix_A9fp10 = (TH2F*)ifile->Get("oP_fpix_regionA_9fp10_METanalysis_Eta2p4");
    TH2F *oP_vs_fpix_D3fp8 = (TH2F*)ifile->Get("oP_fpix_regionD_3fp8_METanalysis_Eta2p4");
    TH2F *oP_vs_fpix_D8fp9 = (TH2F*)ifile->Get("oP_fpix_regionD_8fp9_METanalysis_Eta2p4");
    TH2F *oP_vs_fpix_D9fp10 = (TH2F*)ifile->Get("oP_fpix_regionD_9fp10_METanalysis_Eta2p4");

    TH2D *fpix_vs_ih_C3fp9 = (TH2D*)ifile->Get("fpix_ih_regionC_3fp9_METanalysis_Eta2p4");
    TH2D *fpix_vs_ih_D9fp10_bis = (TH2D*)ifile->Get("fpix_ih_regionD_9fp10_METanalysis_Eta2p4");

    // add the histograms together
    TH2F *ih_vs_eta = (TH2F*)ih_vs_eta_A3fp9->Clone("ih_vs_eta");
    ih_vs_eta->Add(ih_vs_eta_A9fp10);
    ih_vs_eta->Add(ih_vs_eta_D3fp8);
    ih_vs_eta->Add(ih_vs_eta_D8fp9);
    ih_vs_eta->Add(ih_vs_eta_D9fp10);
    ih_vs_eta->RebinX(4);

    TH2F *fpix_vs_ih = (TH2F*)fpix_vs_ih_A3fp9->Clone("fpix_vs_ih");
    fpix_vs_ih->Add(fpix_vs_ih_A9fp10);
    fpix_vs_ih->Add(fpix_vs_ih_D3fp8);
    fpix_vs_ih->Add(fpix_vs_ih_D8fp9);
    fpix_vs_ih->Add(fpix_vs_ih_D9fp10);

    TH2F *oP_vs_fpix = (TH2F*)oP_vs_fpix_A3fp9->Clone("fpix_vs_oP");
    oP_vs_fpix->Add(oP_vs_fpix_A9fp10);
    oP_vs_fpix->Add(oP_vs_fpix_D3fp8);
    oP_vs_fpix->Add(oP_vs_fpix_D8fp9);
    oP_vs_fpix->Add(oP_vs_fpix_D9fp10);
    TH2F *fpix_vs_oP = TransposeTH2(oP_vs_fpix);

    // projection of fpix_vs_ih_C3fp9
    TH1D *py_C3fp9  = fpix_vs_ih_C3fp9->ProjectionY("py_C3fp9");
    TH1D *py_D9fp10 = fpix_vs_ih_D9fp10_bis->ProjectionY("py_D9fp10");
    
    ofile->cd();
    ih_vs_eta->Write();
    fpix_vs_ih->Write();
    fpix_vs_oP->Write();
    pT_vs_fpix->Write();
    py_C3fp9->Write();
    py_D9fp10->Write();
    ofile->Close();

    return;
}

void BKGdependency () {

    TFile *ifile_TTbar = new TFile("../outputDisplay/BKGdependency_preliminary_TTbar.root", "READ");
    TFile *ifile_QCD = new TFile("../outputDisplay/BKGdependency_preliminary_QCD.root", "READ");
    TFile *ifile_Wjets = new TFile("../outputDisplay/BKGdependency_preliminary_Wjets.root", "READ");

    TH2F *ih_vs_eta__TTbar = (TH2F*)ifile_TTbar->Get("ih_vs_eta");
    TH2F *fpix_vs_ih__TTbar = (TH2F*)ifile_TTbar->Get("fpix_vs_ih");
    TH2F *fpix_vs_oP__TTbar = (TH2F*)ifile_TTbar->Get("fpix_vs_oP_transposed");
    TH2F *pT_vs_fpix__TTbar = (TH2F*)ifile_TTbar->Get("METanalysis_Eta2p4_pT_vs_Fpixel");

    TH2F *ih_vs_eta__QCD = (TH2F*)ifile_QCD->Get("ih_vs_eta");
    TH2F *fpix_vs_oP__QCD = (TH2F*)ifile_QCD->Get("fpix_vs_oP_transposed");
    TH2F *fpix_vs_ih__QCD = (TH2F*)ifile_QCD->Get("fpix_vs_ih");
    TH2F *pT_vs_fpix__QCD = (TH2F*)ifile_QCD->Get("METanalysis_Eta2p4_pT_vs_Fpixel");

    TH2F *ih_vs_eta__Wjets = (TH2F*)ifile_Wjets->Get("ih_vs_eta");
    TH2F *fpix_vs_oP__Wjets = (TH2F*)ifile_Wjets->Get("fpix_vs_oP_transposed");
    TH2F *fpix_vs_ih__Wjets = (TH2F*)ifile_Wjets->Get("fpix_vs_ih");
    TH2F *pT_vs_fpix__Wjets = (TH2F*)ifile_Wjets->Get("METanalysis_Eta2p4_pT_vs_Fpixel");

    TH1F *py_C3fp9__TTbar = (TH1F*)ifile_TTbar->Get("py_C3fp9");
    TH1F *py_C3fp9__QCD = (TH1F*)ifile_QCD->Get("py_C3fp9");
    TH1F *py_C3fp9__Wjets = (TH1F*)ifile_Wjets->Get("py_C3fp9");

    TH1F *py_D9fp10__TTbar = (TH1F*)ifile_TTbar->Get("py_D9fp10");
    TH1F *py_D9fp10__QCD = (TH1F*)ifile_QCD->Get("py_D9fp10");
    TH1F *py_D9fp10__Wjets = (TH1F*)ifile_Wjets->Get("py_D9fp10");
    

    // sum the histograms together
    TH2F *ih_vs_eta = (TH2F*)ih_vs_eta__TTbar->Clone("ih_vs_eta_b");
    ih_vs_eta->Add(ih_vs_eta__QCD); ih_vs_eta->Add(ih_vs_eta__Wjets);
    TH2F *fpix_vs_ih = (TH2F*)fpix_vs_ih__TTbar->Clone("fpix_vs_ih_b");
    fpix_vs_ih->Add(fpix_vs_ih__QCD); fpix_vs_ih->Add(fpix_vs_ih__Wjets);
    TH2F *fpix_vs_oP = (TH2F*)fpix_vs_oP__TTbar->Clone("fpix_vs_oP_b");
    fpix_vs_oP->Add(fpix_vs_oP__QCD); fpix_vs_oP->Add(fpix_vs_oP__Wjets);
    TH2F *pT_vs_fpix = (TH2F*)pT_vs_fpix__TTbar->Clone("pT_vs_fpix_b");
    pT_vs_fpix->Add(pT_vs_fpix__QCD); pT_vs_fpix->Add(pT_vs_fpix__Wjets);
    TH1F *py_C3fp9 = (TH1F*)py_C3fp9__TTbar->Clone("py_C3fp9");
    py_C3fp9->Add(py_C3fp9__QCD); py_C3fp9->Add(py_C3fp9__Wjets);
    TH1F *py_D9fp10 = (TH1F*)py_D9fp10__TTbar->Clone("py_D9fp10");
    py_D9fp10->Add(py_D9fp10__QCD); py_D9fp10->Add(py_D9fp10__Wjets);

    // plot the profile
    TProfile *profile_ih_vs_eta = ih_vs_eta->ProfileX("profile_ih_vs_eta");
    TProfile *profile_fpix_vs_ih = fpix_vs_ih->ProfileX("profile_fpix_vs_ih");
    TProfile *profile_fpix_vs_oP = fpix_vs_oP->ProfileX("profile_fpix_vs_oP");
    TProfile *profile_pT_vs_fpix = pT_vs_fpix->ProfileX("profile_pT_vs_fpix");


    fpix_vs_ih->GetXaxis()->SetTitle("F_{pixel}");
    fpix_vs_ih->GetYaxis()->SetTitle("I_{h} [MeV/cm]");
    fpix_vs_ih->GetYaxis()->SetTitleOffset(1.2);
    fpix_vs_ih->GetZaxis()->SetTitle("Events");
    fpix_vs_ih->SetTitle("");
    fpix_vs_ih->GetYaxis()->SetTitleSize(0.06);
    fpix_vs_ih->GetXaxis()->SetTitleSize(0.06);
    fpix_vs_ih->GetXaxis()->SetTitleOffset(0.9);
    fpix_vs_ih->GetYaxis()->SetTitleOffset(1.1);
    fpix_vs_ih->GetXaxis()->SetLabelSize(0.05);
    fpix_vs_ih->GetYaxis()->SetLabelSize(0.05);
    fpix_vs_ih->GetZaxis()->SetTitleSize(0.06);
    fpix_vs_ih->GetZaxis()->SetTitleOffset(0.9);
    fpix_vs_ih->GetZaxis()->SetLabelSize(0.05);

    fpix_vs_oP->GetXaxis()->SetTitle("F_{pixel}");
    fpix_vs_oP->GetYaxis()->SetTitle("10^{4}/p [GeV^{-1}]");
    fpix_vs_oP->GetYaxis()->SetTitleOffset(1.2);
    fpix_vs_oP->GetZaxis()->SetTitle("Events");
    fpix_vs_oP->SetTitle("");
    fpix_vs_oP->GetYaxis()->SetTitleSize(0.06);
    fpix_vs_oP->GetXaxis()->SetTitleSize(0.06);
    fpix_vs_oP->GetXaxis()->SetTitleOffset(0.9);
    fpix_vs_oP->GetYaxis()->SetTitleOffset(1);
    fpix_vs_oP->GetXaxis()->SetLabelSize(0.05);
    fpix_vs_oP->GetYaxis()->SetLabelSize(0.05);
    fpix_vs_oP->GetZaxis()->SetTitleSize(0.06);
    fpix_vs_oP->GetZaxis()->SetTitleOffset(0.9);
    fpix_vs_oP->GetZaxis()->SetLabelSize(0.05);

    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.06);


    TCanvas *c_fpix_vs_ih = new TCanvas("c_fpix_vs_ih", "c_fpix_vs_ih", 800, 600);
    c_fpix_vs_ih->SetLeftMargin(0.16); c_fpix_vs_ih->SetBottomMargin(0.16); c_fpix_vs_ih->SetRightMargin(0.16);
    fpix_vs_ih->Draw("COLZ");
    profile_fpix_vs_ih->SetMarkerStyle(20);
    profile_fpix_vs_ih->SetMarkerColor(kRed);
    profile_fpix_vs_ih->SetLineColor(kRed);
    profile_fpix_vs_ih->Draw("sameP");
    c_fpix_vs_ih->SetLogz();
    gStyle->SetOptStat(0);
    latex1->Draw();
    c_fpix_vs_ih->SaveAs("../outputDisplay/Correlation_fpix_ih.pdf");
    latex1->SetTitle("#it{Private work (CMS simulation)}");
    c_fpix_vs_ih->Modified();
    c_fpix_vs_ih->Update();
    c_fpix_vs_ih->SaveAs("../outputDisplay/Correlation_fpix_ih_bis.pdf");
    latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    fpix_vs_ih->GetYaxis()->SetRangeUser(3.29, 3.31);
    c_fpix_vs_ih->Update();
    c_fpix_vs_ih->SaveAs("../outputDisplay/Correlation_fpix_ih__ZOOM.pdf");
    latex1->SetTitle("#it{Private work (CMS simulation)}");
    c_fpix_vs_ih->Modified();
    c_fpix_vs_ih->Update();
    c_fpix_vs_ih->SaveAs("../outputDisplay/Correlation_fpix_ih__ZOOM_bis.pdf");
    latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");

    TCanvas *c_fpix_vs_oP = new TCanvas("c_fpix_vs_oP", "c_fpix_vs_oP", 800, 600);
    c_fpix_vs_oP->SetLeftMargin(0.16); c_fpix_vs_oP->SetBottomMargin(0.16); c_fpix_vs_oP->SetRightMargin(0.16);
    fpix_vs_oP->Draw("COLZ");
    profile_fpix_vs_oP->SetMarkerStyle(20);
    profile_fpix_vs_oP->SetMarkerColor(kRed);
    profile_fpix_vs_oP->SetLineColor(kRed);
    profile_fpix_vs_oP->Draw("sameP");
    c_fpix_vs_oP->SetLogz();
    gStyle->SetOptStat(0);
    latex1->Draw();
    c_fpix_vs_oP->SaveAs("../outputDisplay/Correlation_fpix_1oP.pdf");
    latex1->SetTitle("#it{Private work (CMS simulation)}");
    c_fpix_vs_oP->Modified();
    c_fpix_vs_oP->Update();
    c_fpix_vs_oP->SaveAs("../outputDisplay/Correlation_fpix_1oP_bis.pdf");
    latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    fpix_vs_oP->GetYaxis()->SetRangeUser(40, 60);
    c_fpix_vs_oP->Update();
    c_fpix_vs_oP->SaveAs("../outputDisplay/Correlation_fpix_1oP__ZOOM.pdf");
    latex1->SetTitle("#it{Private work (CMS simulation)}");
    c_fpix_vs_oP->Modified();
    c_fpix_vs_oP->Update();
    c_fpix_vs_oP->SaveAs("../outputDisplay/Correlation_fpix_1oP__ZOOM_bis.pdf");


    // Ih ratio in C and D
    py_C3fp9->Rebin(4); py_D9fp10->Rebin(4);
    py_C3fp9->Scale(1.0/py_C3fp9->Integral());
    py_D9fp10->Scale(1.0/py_D9fp10->Integral());
    
    py_C3fp9->GetYaxis()->SetTitle("Events (normalised)");
    py_C3fp9->GetXaxis()->SetTitle("I_{h} [MeV/cm]");
    py_C3fp9->SetTitle("");
    py_C3fp9->GetYaxis()->SetTitleSize(0.06);
    py_C3fp9->GetXaxis()->SetTitleSize(0.06);
    py_C3fp9->GetXaxis()->SetTitleOffset(0.9);
    py_C3fp9->GetYaxis()->SetTitleOffset(0.9);
    py_C3fp9->GetXaxis()->SetLabelSize(0.05);
    py_C3fp9->GetYaxis()->SetLabelSize(0.05);
    py_C3fp9->GetXaxis()->SetRangeUser(0, 8);
    py_C3fp9->SetLineColor(kBlue);
    py_C3fp9->SetLineWidth(2);

    py_D9fp10->SetLineColor(kRed);
    py_D9fp10->SetLineWidth(2);
    latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
    pad1->SetLeftMargin(0.16); pad1->SetBottomMargin(0.12);
    pad1->Draw();
    pad1->cd();

    py_C3fp9->Draw("hist");
    py_D9fp10->Draw("hist same");
    pad1->SetLogy();
    
    TLegend *legend = new TLegend(0.5, 0.7, 0.8, 0.9);
    legend->AddEntry(py_C3fp9, "Control region (0.3< F_{pixel}#leq0.8)", "l");
    legend->AddEntry(py_D9fp10, "Validation region (0.8< F_{pixel}#leq0.9)", "l");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.045);
    legend->Draw();
    gStyle->SetOptStat(0);
    latex1->Draw();
    
    // Ratio plot
    c1->cd();
    TPad* pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
    pad2->SetLeftMargin(0.16);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.33);
    pad2->SetTickx(1);
    pad2->SetTicky(1);
    pad2->Draw();
    pad2->cd();

    TH1D *hRatio = (TH1D*)py_C3fp9->Clone("C_over_D");
    hRatio->Divide(py_D9fp10);
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetLineColor(kBlack);


    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->GetXaxis()->SetTitleSize(0.15);
    hRatio->GetXaxis()->SetTitleOffset(0.9);
    hRatio->GetXaxis()->SetLabelSize(0.12);
    hRatio->GetYaxis()->SetTitleSize(0.14);
    hRatio->GetYaxis()->SetTitleOffset(0.3);
    hRatio->GetYaxis()->SetLabelSize(0.12);

    hRatio->Draw("E1");
    TLine *line = new TLine(hRatio->GetXaxis()->GetXmin(), 1, 8, 1);
    line->SetLineStyle(2);
    gStyle->SetOptStat(0);
    line->Draw("same");
    hRatio->GetYaxis()->SetTitle("I_{h}^{CR} / I_{h}^{VR} ");
    hRatio->GetXaxis()->SetTitle("I_{h} [MeV/cm]");
    hRatio->SetMinimum(0.8);
    hRatio->SetMaximum(1.2);

    //fit hRatio with a linear function:
    hRatio->Fit("pol1", "R", "", 2.5, 5);

    c1->SaveAs("../outputDisplay/IhRatio.pdf");
    latex1->SetTitle("#it{Private work (CMS simulation)}");
    c1->Modified();
    c1->Update();
    c1->SaveAs("../outputDisplay/IhRatio_bis.pdf");

    return;
}

void GluinoP_mass(bool isPythia, bool isWeighted) {

    TFile *ofile = new TFile(Form("../outputDisplay/GluinoP_mass_%s_%s.root", isPythia ? "pythia" : "madgraph", isWeighted ? "weighted" : ""), "RECREATE");

    // create a std::vector of ifile where to retrieve the pT vs Fpixel histograms in each:
    std::vector<TFile*> ifiles;
    std::vector<TString> labels;
    if (isPythia) labels = {"1000", "1200", "1400", "1600", "1800", "2000", "2200", "2400", "2600"};
    else labels = {"1100", "1200", "1300", "1400", "1600", "1800", "2000", "2200", "2400", "2600"};

    for (unsigned int i=0; i<labels.size(); i++) ifiles.push_back(new TFile(Form("../output/Gluino_V19/Gluino_Run3_MET_%s%s_V19p%s_%s.root", isPythia ? "" : "madgraph_", labels[i].Data(), isPythia ? "0" : "1", isWeighted ? "weighted" : ""), "READ"));

    // retrieve the p histogram and plot on the same canvas with different color (+legend)
    TCanvas *c_p = new TCanvas("c_p", "c_p", 800, 600);
    TLegend *legend = new TLegend(0.7, 0.5, 0.85, 0.9);
    for (size_t i = 0; i < ifiles.size(); i++) {
        TH1F *p = (TH1F*)ifiles[i]->Get("METanalysis_Eta2p4_P");
        p->SetLineColor(kOrange+i);
        p->SetLineWidth(2);
        if (i == 0) p->Draw("HIST");
        else p->Draw("HIST same");
        legend->AddEntry(p, Form("#tilde{g} m=%s",labels[i].Data()), "l");
    }
    legend->Draw();
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    c_p->SetLogy();

    // retrieve the pT vs Fpixel histograms, sum all the TH2 together and plot it
    TH2F *pT_vs_Fpixel_sum = nullptr;
    float fpix_cut = 0.7;
    std::vector <TH1D*> pT_distributions;
    std::vector <TH1D*> pT_distributions_Fpix_cut;
    std::vector <TH1D*> Fpixel_distributions;
    for (size_t i = 0; i < ifiles.size(); i++) {
        TH2F *pT_vs_Fpixel = (TH2F*)ifiles[i]->Get("METanalysis_Eta2p4_pT_vs_Fpixel");

        float nEntries_CR = 0;
        for (int xbin = 1; xbin <= pT_vs_Fpixel->GetNbinsX(); xbin++) {
            for (int ybin = 1; ybin <= pT_vs_Fpixel->GetNbinsY(); ybin++) {
                if (pT_vs_Fpixel->GetXaxis()->GetBinCenter(xbin) <= fpix_cut || pT_vs_Fpixel->GetYaxis()->GetBinCenter(ybin) <= 70) {
                    nEntries_CR += pT_vs_Fpixel->GetBinContent(xbin, ybin);
                }
            }
        }
        cout << "Number of entries in the CR (pT<=70 and Fpix<=" << fpix_cut << ") for m=" << labels[i].Data() << ": " << nEntries_CR << " / " << pT_vs_Fpixel->Integral() << " (" << 100.0*nEntries_CR/pT_vs_Fpixel->Integral() << "%)" << endl;

        TH1D *proj_pT = pT_vs_Fpixel->ProjectionY(Form("pT_%s", labels[i].Data()));
        TH1D *proj_Fpixel = pT_vs_Fpixel->ProjectionX(Form("Fpixel_%s", labels[i].Data()));
        pT_distributions.push_back(proj_pT);
        Fpixel_distributions.push_back(proj_Fpixel);

        // do the proj of pT but with a Fpix cut on the TH2 before
        TH1D *proj_pT_Fpix_cut = new TH1D(Form("pT_%s_Fpix_cut", labels[i].Data()), Form("pT_%s_Fpix_cut", labels[i].Data()), pT_vs_Fpixel->GetYaxis()->GetNbins(), pT_vs_Fpixel->GetYaxis()->GetXmin(), pT_vs_Fpixel->GetYaxis()->GetXmax());
        for (int xbin = 1; xbin <= pT_vs_Fpixel->GetNbinsX(); xbin++) {
            for (int ybin = 1; ybin <= pT_vs_Fpixel->GetNbinsY(); ybin++) {
                if (pT_vs_Fpixel->GetXaxis()->GetBinCenter(xbin) <= fpix_cut) {
                    proj_pT_Fpix_cut->Fill(pT_vs_Fpixel->GetYaxis()->GetBinCenter(ybin), pT_vs_Fpixel->GetBinContent(xbin, ybin));
                }
            }
        }
        pT_distributions_Fpix_cut.push_back(proj_pT_Fpix_cut);


        if (pT_vs_Fpixel_sum == nullptr) pT_vs_Fpixel_sum = (TH2F*)pT_vs_Fpixel->Clone("pT_vs_Fpixel_sum");
        else pT_vs_Fpixel_sum->Add(pT_vs_Fpixel);
    }
    TCanvas *c_pT_vs_Fpixel = new TCanvas("c_pT_vs_Fpixel", "c_pT_vs_Fpixel", 800, 600);
    pT_vs_Fpixel_sum->GetXaxis()->SetTitle("F_{pixel}");
    pT_vs_Fpixel_sum->GetYaxis()->SetTitle("p_{T} [GeV]");
    pT_vs_Fpixel_sum->Draw("COLZ");
    c_pT_vs_Fpixel->SetLogz();

    float nEntries_CR = 0;
    for (int xbin = 1; xbin <= pT_vs_Fpixel_sum->GetNbinsX(); xbin++) {
        for (int ybin = 1; ybin <= pT_vs_Fpixel_sum->GetNbinsY(); ybin++) {
            if (pT_vs_Fpixel_sum->GetXaxis()->GetBinCenter(xbin) <= fpix_cut || pT_vs_Fpixel_sum->GetYaxis()->GetBinCenter(ybin) <= 70) {
                nEntries_CR += pT_vs_Fpixel_sum->GetBinContent(xbin, ybin);
            }
        }
    }
    cout << "Number of entries in the CR (pT<=70 and Fpix<=" << fpix_cut << ") (TOTAL): " << nEntries_CR << " / " << pT_vs_Fpixel_sum->Integral() << " (" << 100.0*nEntries_CR/pT_vs_Fpixel_sum->Integral() << "%)" << endl;

    // plot pT_distributions on the same canvas with different color (+legend)
    TCanvas *c_pT = new TCanvas("c_pT", "c_pT", 800, 600);
    TLegend *legend_pT = new TLegend(0.7, 0.5, 0.85, 0.9);
    for (size_t i = 0; i < pT_distributions.size(); i++) {
        pT_distributions[i]->SetLineColor(kOrange+i);
        pT_distributions[i]->SetLineWidth(2);
        if (i == 0) pT_distributions[i]->Draw("HIST");
        else pT_distributions[i]->Draw("HIST same");
        legend_pT->AddEntry(pT_distributions[i], Form("#tilde{g} m=%s",labels[i].Data()), "l");
    }
    legend_pT->Draw();
    legend_pT->SetBorderSize(0);
    legend_pT->SetFillStyle(0);

    // plot Fpixel_distributions on the same canvas with different color (+legend)
    TCanvas *c_Fpixel = new TCanvas("c_Fpixel", "c_Fpixel", 800, 600);
    TLegend *legend_Fpixel = new TLegend(0.7, 0.5, 0.85, 0.9);
    for (size_t i = 0; i < Fpixel_distributions.size(); i++) {
        Fpixel_distributions[i]->SetLineColor(kOrange+i);
        Fpixel_distributions[i]->SetLineWidth(2);
        if (i == 0) Fpixel_distributions[i]->Draw("HIST");
        else Fpixel_distributions[i]->Draw("HIST same");
        legend_Fpixel->AddEntry(Fpixel_distributions[i], Form("#tilde{g} m=%s",labels[i].Data()), "l");
    }
    legend_Fpixel->Draw();
    legend_Fpixel->SetBorderSize(0);
    legend_Fpixel->SetFillStyle(0);

    TCanvas *c_pT2 = new TCanvas("c_pT2", "c_pT2", 800, 600);
    for (size_t i = 0; i < pT_distributions_Fpix_cut.size(); i++) {
        pT_distributions_Fpix_cut[i]->SetLineColor(kOrange+i);
        pT_distributions_Fpix_cut[i]->SetLineWidth(2);
        if (i == 0) pT_distributions_Fpix_cut[i]->Draw("HIST");
        else pT_distributions_Fpix_cut[i]->Draw("HIST same");
    }
    legend_pT->Draw();



    ofile->cd();
    c_p->Write();
    pT_vs_Fpixel_sum->Write();
    c_pT->Write();
    c_Fpixel->Write();
    c_pT2->Write();
    ofile->Close();

}

void CompareMaping(const char *ttbar = "../output/TTbar2024_V15/TTbar2024_V15p9_weighted.root",
                   const char *wjets = "../output/Wjets2024_V14/WjetMuNu2024_V14p13_weighted.root",
                   const char *qcd = "../output/QCD2024_V16/QCD2024_mu_V16p2_weighted.root",
                   const char *gluino = "../output/Gluino_V19/Gluino_Run3_MET_madgraph_2000_V19p8_weighted.root") {

    TFile *ofile = new TFile(OutPath("CompareMapping.root"), "RECREATE");

    // backgrounds: an empty path leaves the sample out of the sum
    TFile *ifileTTbar2024 = OpenOptional(ttbar);
    TFile *ifileWjets2024 = OpenOptional(wjets);
    TFile *ifileQCD2024 = OpenOptional(qcd);
    TFile *GluinoRun3madgraph = new TFile(gluino, "READ");

    TLatex *latex1 = new TLatex(0.155, 0.91, "#it{Private work (CMS simulation)}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);


    // load histograms
    TH2F *pT_vs_fpix_gluino = (TH2F*)GluinoRun3madgraph->Get("METanalysis_TestPUppiMETCut_Eta2p4_pT_vs_Fpixel");

    // sum the MC backgrounds given together
    TH2F *pT_vs_fpix_bkg = nullptr;
    for (TFile *f : {ifileTTbar2024, ifileWjets2024, ifileQCD2024}) {
        if (!f) continue;
        TH2F *h = (TH2F*)f->Get("METanalysis_TestPUppiMETCut_Eta2p4_pT_vs_Fpixel");
        if (!h) Abort(Form("CompareMaping: no METanalysis_TestPUppiMETCut_Eta2p4_pT_vs_Fpixel in %s", f->GetName()));
        if (!pT_vs_fpix_bkg) pT_vs_fpix_bkg = (TH2F*)h->Clone("pT_vs_fpix_bkg");
        else pT_vs_fpix_bkg->Add(h);
    }
    if (!pT_vs_fpix_bkg) Abort("CompareMaping: give at least one background");

    // divide the gluino histogram by the bkg histogram to see where the signal is enhanced
    TH2F *pT_vs_fpix_ratio = (TH2F*)pT_vs_fpix_gluino->Clone("pT_vs_fpix_ratio");
    pT_vs_fpix_ratio->Divide(pT_vs_fpix_bkg);
    pT_vs_fpix_ratio->Rebin2D(5,5);

    // canvas of the ratio
    TCanvas *c_ratio = new TCanvas("c_ratio", "c_ratio", 800, 600);
    pT_vs_fpix_ratio->GetXaxis()->SetTitle("F_{pixel}");
    pT_vs_fpix_ratio->GetYaxis()->SetTitle("p_{T} [GeV]");
    pT_vs_fpix_ratio->Draw("TEXT");

    TCanvas *c_bkg_gluino = new TCanvas("c_bkg_gluino", "c_bkg_gluino", 800, 600);
    c_bkg_gluino->SetLeftMargin(0.16); c_bkg_gluino->SetRightMargin(0.16); c_bkg_gluino->SetBottomMargin(0.16);
    pT_vs_fpix_bkg->GetXaxis()->SetTitle("F_{pixel}");
    pT_vs_fpix_bkg->GetYaxis()->SetTitle("p_{T} [GeV]");
    pT_vs_fpix_bkg->Draw("COLZ");

    TProfile *profile_pT_vs_fpix_bkg = pT_vs_fpix_bkg->ProfileX("profile_pT_vs_fpix_bkg");
    profile_pT_vs_fpix_bkg->SetLineColor(kRed);
    profile_pT_vs_fpix_bkg->SetMarkerColor(kRed);
    profile_pT_vs_fpix_bkg->Draw("sameP");

    // save the profile in a .root file:
    

    pT_vs_fpix_bkg->GetYaxis()->SetTitleSize(0.06);
    pT_vs_fpix_bkg->GetXaxis()->SetTitleSize(0.06);
    pT_vs_fpix_bkg->GetXaxis()->SetTitleOffset(0.9);
    pT_vs_fpix_bkg->GetYaxis()->SetTitleOffset(1.1);
    pT_vs_fpix_bkg->GetXaxis()->SetLabelSize(0.05);
    pT_vs_fpix_bkg->GetYaxis()->SetLabelSize(0.05);
    pT_vs_fpix_bkg->GetZaxis()->SetTitleSize(0.06);
    pT_vs_fpix_bkg->GetZaxis()->SetTitleOffset(0.9);
    pT_vs_fpix_bkg->GetZaxis()->SetLabelSize(0.05);

    pT_vs_fpix_gluino->GetYaxis()->SetTitleSize(0.06);
    pT_vs_fpix_gluino->GetXaxis()->SetTitleSize(0.06);
    pT_vs_fpix_gluino->GetXaxis()->SetTitleOffset(0.9);
    pT_vs_fpix_gluino->GetYaxis()->SetTitleOffset(1.1);
    pT_vs_fpix_gluino->GetXaxis()->SetLabelSize(0.05);
    pT_vs_fpix_gluino->GetYaxis()->SetLabelSize(0.05);
    pT_vs_fpix_gluino->GetZaxis()->SetTitleSize(0.06);
    pT_vs_fpix_gluino->GetZaxis()->SetTitleOffset(0.9);
    pT_vs_fpix_gluino->GetZaxis()->SetLabelSize(0.05);

    pT_vs_fpix_gluino->GetXaxis()->SetTitle("F_{pixel}");
    pT_vs_fpix_gluino->GetYaxis()->SetTitle("p_{T} [GeV]");
    pT_vs_fpix_gluino->GetZaxis()->SetTitle("Entries");

    pT_vs_fpix_bkg->GetXaxis()->SetTitle("F_{pixel}");
    pT_vs_fpix_bkg->GetYaxis()->SetTitle("p_{T} [GeV]");
    pT_vs_fpix_bkg->GetZaxis()->SetTitle("Entries");

    gStyle->SetPalette(kViridis);

    //pT_vs_fpix_gluino->Draw("COLZ");
    pT_vs_fpix_bkg->Draw("COLZ");

    profile_pT_vs_fpix_bkg->Draw("sameP");

    profile_pT_vs_fpix_bkg->Fit("pol1", "R", "", 0.3, 1);
    double r = pT_vs_fpix_bkg->GetCorrelationFactor();
    std::cout << "Correlation factor pT vs Fpix: " << r << std::endl;

    latex1->Draw("same");
    pT_vs_fpix_bkg->GetYaxis()->SetRangeUser(100,150);
    c_bkg_gluino->SetLogz();

    //TLatex *mtext = new TLatex(0.68, 0.91, "#scale[1.3]{#bf{m_{#tilde{g}}=2000 GeV}}");
    TLatex *mtext = new TLatex(0.72, 0.91, "#scale[1.3]{#bf{SM MC}}");
    mtext->SetNDC(); mtext->SetTextFont(42); mtext->SetTextSize(0.04);
    mtext->Draw();

    const double pTcut = 70.;
    const int nY  = pT_vs_fpix_bkg->GetNbinsY();
    const int yLo = pT_vs_fpix_bkg->GetYaxis()->FindBin(pTcut);
    std::cout << "Bord bas du bin de coupure : "
            << pT_vs_fpix_bkg->GetYaxis()->GetBinLowEdge(yLo) << std::endl;  // doit valoir 70

    TH1D *hAll  = pT_vs_fpix_bkg->ProjectionX("hAll",  1,   nY+1);  // overflow inclus
    TH1D *hPass = pT_vs_fpix_bkg->ProjectionX("hPass", yLo, nY+1);

    //hAll->Rebin(4); hPass->Rebin(4);   // à ajuster selon la statistique

    TH1D *hEff = (TH1D*)hPass->Clone("hEff");
    hEff->Divide(hPass, hAll, 1., 1., "B");
    hEff->GetXaxis()->SetTitle("F_{pixel}");
    hEff->GetYaxis()->SetTitle(Form("#varepsilon(p_{T} > %.0f GeV)", pTcut));

    TF1 *f0 = new TF1("f0", "pol0", 0., 1.);
    hEff->Fit(f0, "R");
    std::cout << "chi2/ndf = " << f0->GetChisquare() / f0->GetNDF()
            << "   p-value = " << f0->GetProb() << std::endl;

    // 1. Borne sur une dérive linéaire résiduelle
    TF1 *f1 = new TF1("f1", "pol1", 0., 1.);
    hEff->Fit(f1, "R+");
    std::cout << "pente = " << f1->GetParameter(1) << " +/- " << f1->GetParError(1) << std::endl;

    // 2. Double rapport à la frontière ABCD (Fpixel = 0.9)
    TH1D *hFail = pT_vs_fpix_bkg->ProjectionX("hFail", 1, yLo-1);

    auto effErr = [&](int b1, int b2, double &eps, double &sig) {
        double sP, sF;
        double P = hPass->IntegralAndError(b1, b2, sP);
        double F = hFail->IntegralAndError(b1, b2, sF);
        double A = P + F;
        eps = P / A;
        sig = sqrt((1-eps)*(1-eps)*sP*sP + eps*eps*sF*sF) / A;
    };

    int xB = hEff->GetXaxis()->FindBin(0.9);
    double eLo, sLo, eHi, sHi;
    effErr(1, xB-1, eLo, sLo);
    effErr(xB, hAll->GetNbinsX()+1, eHi, sHi);

    double D  = eHi/eLo - 1.;
    double sD = (eHi/eLo) * sqrt(pow(sHi/eHi,2) + pow(sLo/eLo,2));
    std::cout << "Delta = " << 100*D << " +/- " << 100*sD << " %" << std::endl;


    // cout number of event bellow each .1 in Fpixel, for both hists.
    
    for (int i=1; i<=10; i++) {
        float fpix_cut = 0.1*i;
        int nEntries_bkg = 0;
        float nEntries_gluino = 0;
        for (int xbin = 1; xbin <= pT_vs_fpix_bkg->GetNbinsX(); xbin++) {
            for (int ybin = 1; ybin <= pT_vs_fpix_bkg->GetNbinsY(); ybin++) {
                if (pT_vs_fpix_bkg->GetXaxis()->GetBinCenter(xbin) <= fpix_cut) {
                    nEntries_bkg += pT_vs_fpix_bkg->GetBinContent(xbin, ybin);
                    nEntries_gluino += pT_vs_fpix_gluino->GetBinContent(xbin, ybin);
                }
            }
        }
        cout << "Number of entries with Fpixel <= " << fpix_cut << ": " << nEntries_bkg << " (bkg), " << nEntries_gluino << " (gluino)" << "    (" << 100.0*nEntries_gluino/nEntries_bkg << " %)" << endl;
    }
    gStyle->SetOptStat(0);
    c_bkg_gluino->SaveAs(OutPath("CompareMaping_prout.pdf"));


    ofile->cd();
    c_ratio->Write();
    c_bkg_gluino->Write();
    profile_pT_vs_fpix_bkg->Write();
    ofile->Close();

    return;
}

void FpixelPlot(const char *wjets = "../output/Wjets2024_V14/WjetMuNu2024_V14p11_weighted.root",
                const char *data = "../output/JetMET2024_V12/JetMET2024_V12p24.root",
                const char *gluino = "../output/Gluino_V19/Gluino_Run3_MET_madgraph_2000_V19p1.root") {

    TFile *ofile = new TFile(OutPath("Nm1_plots.root"), "RECREATE");

    TFile *fileWjets = OpenOptional(wjets);                // the background may be left out (empty path)
    TFile *ifile_JetMETdata = new TFile(data, "READ");
    TFile *ifile_Gluino = new TFile(gluino, "READ");


    TH1F *fpix_JetMETdata = (TH1F*)ifile_JetMETdata->Get("Nosel_Fpix");
    TH1F *fpix_Gluino = (TH1F*)ifile_Gluino->Get("Nosel_Fpix");
    TH1F *fpix_Wjets = fileWjets ? (TH1F*)fileWjets->Get("Nosel_Fpix") : nullptr;
    fpix_JetMETdata->Scale(1./fpix_JetMETdata->Integral());
    fpix_Gluino->Scale(1./fpix_Gluino->Integral());
    if (fpix_Wjets) {
        fpix_Wjets->Scale(1./fpix_Wjets->Integral());
        for (int bin = 1; bin <= 5; bin++) fpix_Wjets->SetBinContent(bin, 0);
        fpix_Wjets->SetLineWidth(2);
        fpix_Wjets->SetLineColor(kBlue);
    }

    fpix_Gluino->SetBinContent(1,0);
    fpix_Gluino->SetBinContent(2,0);
    fpix_Gluino->SetBinContent(3,0);
    fpix_Gluino->SetBinContent(4,0);
    fpix_Gluino->SetBinContent(5,0);

    //draw vertical dashed line black at fpix = 0.9
    TLine *line = new TLine(0.9, 0, 0.9, fpix_Gluino->GetMaximum()*1.05);
    line->SetLineStyle(2);
    line->SetLineColor(kBlack);
    line->SetLineWidth(2);
    


    fpix_Gluino->GetYaxis()->SetTitleSize(0.06);
    fpix_Gluino->GetXaxis()->SetTitleSize(0.06);
    fpix_Gluino->GetXaxis()->SetTitleOffset(0.9);
    fpix_Gluino->GetYaxis()->SetTitleOffset(1.1);
    fpix_Gluino->GetXaxis()->SetLabelSize(0.05);
    fpix_Gluino->GetYaxis()->SetLabelSize(0.05);

    TCanvas *c_fpix = new TCanvas("c_fpix", "c_fpix", 800, 600);
    c_fpix->SetLeftMargin(0.16); c_fpix->SetBottomMargin(0.16);
    fpix_JetMETdata->SetLineColor(kBlack);
    fpix_JetMETdata->SetMarkerColor(kBlack);
    fpix_JetMETdata->SetMarkerStyle(20);
    fpix_Gluino->GetXaxis()->SetTitle("F_{pixel}");
    fpix_Gluino->GetYaxis()->SetTitle("Events (normalised)");
    fpix_Gluino->SetFillColorAlpha(kRed, 0.5);
    fpix_Gluino->SetLineColor(kRed);
    fpix_Gluino->Draw("HIST");
    if (fpix_Wjets) fpix_Wjets->Draw("HIST same");
    line->Draw("same");
    //fpix_JetMETdata->Draw("E1 same");
    TLegend *legend_fpix = new TLegend(0.2, 0.5, 0.5, 0.7);
    //legend_fpix->AddEntry(fpix_JetMETdata, "JetMET 2024 data", "lep");
    if (fpix_Wjets) legend_fpix->AddEntry(fpix_Wjets, Label("wjets", "MC Wjet 2024"), "lep");
    legend_fpix->AddEntry(fpix_Gluino, Label("gluino", "#tilde{g}_{m=2000 GeV}"), "f");
    legend_fpix->SetBorderSize(0);
    legend_fpix->SetFillStyle(0);
    legend_fpix->Draw();

    // remove stat using SetOptStat(0);
    gStyle->SetOptStat(0);
    
    TLatex *latex1 = new TLatex(0.155, 0.91, "#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);
    latex1->Draw();
    //latex1->SetTitle("#it{Private work (CMS simulation)}");

    c_fpix->SaveAs(OutPath("Fpixel_plot.pdf"));
}


void PseudoMET_vs_PFMET () {
    TFile *ifileData = new TFile("TestMET2024_V12p24.root", "READ");
    TFile *ifileGluino = new TFile("../output/Gluino_V19/Gluino_Run3_MET_madgraph_2000_V19p1.root", "READ");

    TH2F *CaloMET_vs_PFMET_data = (TH2F*)ifileData->Get("Nosel_PseudoMET_vs_PFMET");
    TH2F *CaloMET_vs_PFMET_gluino = (TH2F*)ifileGluino->Get("Nosel_PseudoMET_vs_PFMET");

    // draw two different TCanvas with the same range, one for data and one for gluino, with the same color scale
    TCanvas *c_data = new TCanvas("c_data", "c_data", 1200, 900);
    CaloMET_vs_PFMET_data->GetXaxis()->SetTitle("PFMET [GeV]");
    CaloMET_vs_PFMET_data->GetXaxis()->SetLabelSize(0.03);
    CaloMET_vs_PFMET_data->GetYaxis()->SetTitle("PseudoMET [GeV]");
    CaloMET_vs_PFMET_data->GetYaxis()->SetLabelSize(0.03);
    CaloMET_vs_PFMET_data->Draw("COLZ");
    c_data->SetLogz();
    gStyle->SetOptStat(0);
    TCanvas *c_gluino = new TCanvas("c_gluino", "c_gluino", 1200, 900);
    CaloMET_vs_PFMET_gluino->GetXaxis()->SetTitle("PFMET [GeV]");
    CaloMET_vs_PFMET_gluino->GetXaxis()->SetLabelSize(0.03);
    CaloMET_vs_PFMET_gluino->GetYaxis()->SetTitle("PseudoMET [GeV]");
    CaloMET_vs_PFMET_gluino->GetYaxis()->SetLabelSize(0.03);
    CaloMET_vs_PFMET_gluino->Draw("COLZ");
    c_gluino->SetLogz();
    gStyle->SetOptStat(0);

    c_data->SaveAs("../outputDisplay/PseudoMET_vs_PFMET_data.pdf");
    c_gluino->SaveAs("../outputDisplay/PseudoMET_vs_PFMET_gluino.pdf");

    return;
}

void DrawPseudoMET(bool OTHERrescale = false) {

    gErrorIgnoreLevel = kError;

    TFile *ifileWjetMuNu;
    if (OTHERrescale) ifileWjetMuNu = new TFile("../output/Wjets2024_V14/WjetMuNu2024_V14p10_weighted.root", "READ");
    else ifileWjetMuNu = new TFile("../output/Wjets2024_V14/WjetMuNu2024_V14p8_weighted.root", "READ");
    TFile *ifileMET = new TFile("../output/JetMET2024_V12/JetMET2024_V12p24.root", "READ");

    TH1F *Nosel_PseudoCaloMET_WjetMuNu = (TH1F*)ifileWjetMuNu->Get("Nosel_PseudoCaloMET");
    TH1F *Nosel_PseudoCaloMET_MET2024 = (TH1F*)ifileMET->Get("Nosel_PseudoCaloMET");

    TH1F *Nm1_event_PseudoMET_WjetMuNu = (TH1F*)ifileWjetMuNu->Get("Nm1_event_CaloMET");
    TH1F *Nm1_event_PseudoMET_MET2024 = (TH1F*)ifileMET->Get("Nm1_event_CaloMET");

    TH1F *Nm1Rescaled_event_PseudoMET_WjetMuNu = (TH1F*)ifileWjetMuNu->Get(Form("Nm1_event_CaloMET_%s", OTHERrescale ? "rescaled" : "weighted"));
    TH1F *Nm1Rescaled_event_PseudoMET_MET2024 = (TH1F*)ifileMET->Get("Nm1_event_CaloMET_weighted");


    Nm1_event_PseudoMET_WjetMuNu->Scale(1./Nm1_event_PseudoMET_WjetMuNu->Integral());
    Nm1_event_PseudoMET_MET2024->Scale(1./Nm1_event_PseudoMET_MET2024->Integral());
    Nosel_PseudoCaloMET_WjetMuNu->Scale(1./Nosel_PseudoCaloMET_WjetMuNu->Integral());
    Nosel_PseudoCaloMET_MET2024->Scale(1./Nosel_PseudoCaloMET_MET2024->Integral());
    Nm1Rescaled_event_PseudoMET_WjetMuNu->Scale(1./Nm1Rescaled_event_PseudoMET_WjetMuNu->Integral());
    Nm1Rescaled_event_PseudoMET_MET2024->Scale(1./Nm1Rescaled_event_PseudoMET_MET2024->Integral());

    Nosel_PseudoCaloMET_WjetMuNu->GetXaxis()->SetTitle("Pseudo MET [GeV]");
    Nosel_PseudoCaloMET_WjetMuNu->GetYaxis()->SetTitle("Number of events");
    Nosel_PseudoCaloMET_WjetMuNu->GetYaxis()->SetTitleOffset(1.2);
    Nosel_PseudoCaloMET_WjetMuNu->SetTitle("");
    Nosel_PseudoCaloMET_WjetMuNu->GetYaxis()->SetTitleSize(0.06);
    Nosel_PseudoCaloMET_WjetMuNu->GetXaxis()->SetTitleSize(0.06);
    Nosel_PseudoCaloMET_WjetMuNu->GetXaxis()->SetTitleOffset(0.9);
    Nosel_PseudoCaloMET_WjetMuNu->GetYaxis()->SetTitleOffset(1);
    Nosel_PseudoCaloMET_WjetMuNu->GetXaxis()->SetLabelSize(0.05);
    Nosel_PseudoCaloMET_WjetMuNu->GetYaxis()->SetLabelSize(0.05);

    Nm1_event_PseudoMET_WjetMuNu->GetXaxis()->SetTitle("Pseudo MET [GeV]");
    Nm1_event_PseudoMET_WjetMuNu->GetYaxis()->SetTitle("Number of events");
    Nm1_event_PseudoMET_WjetMuNu->GetYaxis()->SetTitleOffset(1.2);
    Nm1_event_PseudoMET_WjetMuNu->SetTitle("");
    Nm1_event_PseudoMET_WjetMuNu->GetYaxis()->SetTitleSize(0.06);
    Nm1_event_PseudoMET_WjetMuNu->GetXaxis()->SetTitleSize(0.06);
    Nm1_event_PseudoMET_WjetMuNu->GetXaxis()->SetTitleOffset(0.9);
    Nm1_event_PseudoMET_WjetMuNu->GetYaxis()->SetTitleOffset(1);
    Nm1_event_PseudoMET_WjetMuNu->GetXaxis()->SetLabelSize(0.05);
    Nm1_event_PseudoMET_WjetMuNu->GetYaxis()->SetLabelSize(0.05);

    Nm1Rescaled_event_PseudoMET_WjetMuNu->GetXaxis()->SetTitle("Pseudo MET [GeV]");
    Nm1Rescaled_event_PseudoMET_WjetMuNu->GetYaxis()->SetTitle("Number of events");
    Nm1Rescaled_event_PseudoMET_WjetMuNu->GetYaxis()->SetTitleOffset(1.2);
    Nm1Rescaled_event_PseudoMET_WjetMuNu->SetTitle("");
    Nm1Rescaled_event_PseudoMET_WjetMuNu->GetYaxis()->SetTitleSize(0.06);
    Nm1Rescaled_event_PseudoMET_WjetMuNu->GetXaxis()->SetTitleSize(0.06);
    Nm1Rescaled_event_PseudoMET_WjetMuNu->GetXaxis()->SetTitleOffset(0.9);
    Nm1Rescaled_event_PseudoMET_WjetMuNu->GetYaxis()->SetTitleOffset(1);
    Nm1Rescaled_event_PseudoMET_WjetMuNu->GetXaxis()->SetLabelSize(0.05);
    Nm1Rescaled_event_PseudoMET_WjetMuNu->GetYaxis()->SetLabelSize(0.05);

    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);

    // fit both hist with a gaussian at the peak
    TF1 *gaus_WjetMuNu = new TF1("gaus_WjetMuNu", "gaus", 20, 200);
    TF1 *gaus_MET2024 = new TF1("gaus_MET2024", "gaus", 20, 200);
    gaus_WjetMuNu->SetParameters(Nosel_PseudoCaloMET_WjetMuNu->GetMaximum(), Nosel_PseudoCaloMET_WjetMuNu->GetMean(), Nosel_PseudoCaloMET_WjetMuNu->GetRMS());
    gaus_MET2024->SetParameters(Nosel_PseudoCaloMET_MET2024->GetMaximum(), Nosel_PseudoCaloMET_MET2024->GetMean(), Nosel_PseudoCaloMET_MET2024->GetRMS());
    cout << "No sel" << endl;
    Nosel_PseudoCaloMET_WjetMuNu->Fit(gaus_WjetMuNu, "RQ");
    Nosel_PseudoCaloMET_MET2024->Fit(gaus_MET2024, "RQ");
    gaus_WjetMuNu->SetLineColor(kRed);
    gaus_MET2024->SetLineColor(kGreen);

    TF1 *gaus_Nm1_WjetMuNu = new TF1("gaus_Nm1_WjetMuNu", "gaus", 80, 220);
    TF1 *gaus_Nm1_MET2024 = new TF1("gaus_Nm1_MET2024", "gaus", 50, 230);
    gaus_Nm1_WjetMuNu->SetParameters(Nm1_event_PseudoMET_WjetMuNu->GetMaximum(), Nm1_event_PseudoMET_WjetMuNu->GetMean(), Nm1_event_PseudoMET_WjetMuNu->GetRMS());
    gaus_Nm1_MET2024->SetParameters(Nm1_event_PseudoMET_MET2024->GetMaximum(), Nm1_event_PseudoMET_MET2024->GetMean(), Nm1_event_PseudoMET_MET2024->GetRMS());
    cout << "N-1 sel" << endl;
    Nm1_event_PseudoMET_WjetMuNu->Fit(gaus_Nm1_WjetMuNu, "R");
    Nm1_event_PseudoMET_MET2024->Fit(gaus_Nm1_MET2024, "R");
    gaus_Nm1_WjetMuNu->SetLineColor(kRed);
    gaus_Nm1_MET2024->SetLineColor(kGreen);
    cout << "Mean data: " << gaus_Nm1_MET2024->GetParameter(1) << "   Mean WjetMuNu: " << gaus_Nm1_WjetMuNu->GetParameter(1) << endl;
    cout << "weight = " << gaus_Nm1_MET2024->GetParameter(1)/gaus_Nm1_WjetMuNu->GetParameter(1) << endl;

    TLegend *legend = new TLegend(0.7, 0.5, 0.85, 0.9);
    legend->AddEntry(Nosel_PseudoCaloMET_WjetMuNu, "WjetMuNu", "f");
    legend->AddEntry(Nosel_PseudoCaloMET_MET2024, "MET2024", "lep");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);


    TCanvas *c_PseudoCaloMET_nosel = new TCanvas("c_PseudoCaloMET_nosel", "c_PseudoCaloMET_nosel", 800, 600);
    c_PseudoCaloMET_nosel->cd();
    c_PseudoCaloMET_nosel->SetLeftMargin(0.16);
    c_PseudoCaloMET_nosel->SetBottomMargin(0.16);
    Nosel_PseudoCaloMET_WjetMuNu->SetLineColor(kBlue-7);
    Nosel_PseudoCaloMET_WjetMuNu->SetFillColorAlpha(kBlue-7, 0.5);
    Nosel_PseudoCaloMET_WjetMuNu->Draw("HIST");
    Nosel_PseudoCaloMET_WjetMuNu->GetXaxis()->SetRangeUser(10, 1000);
    Nosel_PseudoCaloMET_MET2024->SetLineColor(kBlack);
    Nosel_PseudoCaloMET_MET2024->SetMarkerColor(kBlack);
    Nosel_PseudoCaloMET_MET2024->SetMarkerStyle(20);
    Nosel_PseudoCaloMET_MET2024->Draw("E1 same");
    gaus_WjetMuNu->Draw("same");
    gaus_MET2024->Draw("same");
    latex1->Draw();
    legend->Draw();
    gPad->SetLogy();
    gStyle->SetOptStat(0);

    TCanvas *c_PseudoCaloMET = new TCanvas("c_PseudoCaloMET", "c_PseudoCaloMET", 800, 600);
    c_PseudoCaloMET->cd();
    c_PseudoCaloMET->SetLeftMargin(0.16);
    c_PseudoCaloMET->SetBottomMargin(0.16);
    Nm1_event_PseudoMET_WjetMuNu->SetLineColor(kBlue-7);
    Nm1_event_PseudoMET_WjetMuNu->SetFillColorAlpha(kBlue-7, 0.5);
    Nm1_event_PseudoMET_WjetMuNu->Draw("HIST");
    Nm1_event_PseudoMET_WjetMuNu->GetXaxis()->SetRangeUser(10, 1000);
    Nm1_event_PseudoMET_MET2024->SetLineColor(kBlack);
    Nm1_event_PseudoMET_MET2024->SetMarkerColor(kBlack);
    Nm1_event_PseudoMET_MET2024->SetMarkerStyle(20);
    Nm1_event_PseudoMET_MET2024->Draw("E1 same");
    gaus_Nm1_WjetMuNu->Draw("same");
    gaus_Nm1_MET2024->Draw("same");
    legend->Draw();
    latex1->Draw();
    gPad->SetLogy();
    gStyle->SetOptStat(0);

    TCanvas *cRescaled_PseudoCaloMET = new TCanvas("cRescaled_PseudoCaloMET", "cRescaled_PseudoCaloMET", 800, 600);
    cRescaled_PseudoCaloMET->cd();
    cRescaled_PseudoCaloMET->SetLeftMargin(0.16);
    cRescaled_PseudoCaloMET->SetBottomMargin(0.16);
    Nm1Rescaled_event_PseudoMET_WjetMuNu->SetLineColor(kBlue-7);
    Nm1Rescaled_event_PseudoMET_WjetMuNu->SetFillColorAlpha(kBlue-7, 0.5);
    Nm1Rescaled_event_PseudoMET_WjetMuNu->Draw("HIST");
    Nm1Rescaled_event_PseudoMET_WjetMuNu->GetXaxis()->SetRangeUser(10, 1000);
    Nm1Rescaled_event_PseudoMET_MET2024->SetLineColor(kBlack);
    Nm1Rescaled_event_PseudoMET_MET2024->SetMarkerColor(kBlack);
    Nm1Rescaled_event_PseudoMET_MET2024->SetMarkerStyle(20);
    Nm1Rescaled_event_PseudoMET_MET2024->Draw("E1 same");
    legend->Draw();
    latex1->Draw();
    gPad->SetLogy();
    gStyle->SetOptStat(0);

    TCanvas *cRatio_nosel = DrawWithRatio(Nosel_PseudoCaloMET_MET2024, Nosel_PseudoCaloMET_WjetMuNu,
                       c_PseudoCaloMET_nosel, "PseudoMET wo selections", "data/MC",  "Pseudo MET [GeV]", 10, 1000);

    TCanvas *cRatio_Nm1 = DrawWithRatio(Nm1_event_PseudoMET_MET2024, Nm1_event_PseudoMET_WjetMuNu,
                       c_PseudoCaloMET, "PseudoMET N-1 selections", "data/MC",  "Pseudo MET [GeV]", 10, 1000);

    TCanvas *cRatio_Rescaled = DrawWithRatio(Nm1Rescaled_event_PseudoMET_MET2024, Nm1Rescaled_event_PseudoMET_WjetMuNu,
                       cRescaled_PseudoCaloMET, "PseudoMET N-1 selections rescaled", "data/MC",  "Pseudo MET [GeV]", 10, 1000);

    cRatio_nosel->SaveAs("../outputDisplay/PseudoCaloMET_nosel_FIT.pdf");
    cRatio_Nm1->SaveAs("../outputDisplay/PseudoCaloMET_Nm1_FIT.pdf");
    cRatio_Rescaled->SaveAs(Form("../outputDisplay/PseudoCaloMET_Nm1_%srescaled.pdf", OTHERrescale ? "OTHER" : ""));


    latex1->SetTitle("#it{Private work (CMS simulation/data)}");
    
    c_PseudoCaloMET_nosel->Modified();
    c_PseudoCaloMET_nosel->Update();
    TCanvas *cRatio_nosel_bis = DrawWithRatio(Nosel_PseudoCaloMET_MET2024, Nosel_PseudoCaloMET_WjetMuNu,
                       c_PseudoCaloMET_nosel, "PseudoMET wo selections", "data/MC",  "Pseudo MET [GeV]", 10, 1000);
    cRatio_nosel_bis->SaveAs("../outputDisplay/IhRatio_bis.pdf");

    c_PseudoCaloMET->Modified();
    c_PseudoCaloMET->Update();
    TCanvas *cRatio_Nm1_bis = DrawWithRatio(Nm1_event_PseudoMET_MET2024, Nm1_event_PseudoMET_WjetMuNu,
                       c_PseudoCaloMET, "PseudoMET N-1 selections", "data/MC",  "Pseudo MET [GeV]", 10, 1000);
    cRatio_Nm1_bis->SaveAs("../outputDisplay/PseudoCaloMET_Nm1_FIT_bis.pdf");

    cRescaled_PseudoCaloMET->Modified();
    cRescaled_PseudoCaloMET->Update();
    TCanvas *cRatio_Rescaled_bis = DrawWithRatio(Nm1Rescaled_event_PseudoMET_MET2024, Nm1Rescaled_event_PseudoMET_WjetMuNu,
                       cRescaled_PseudoCaloMET, "PseudoMET N-1 selections rescaled", "data/MC",  "Pseudo MET [GeV]", 10, 1000);
    cRatio_Rescaled_bis->SaveAs(Form("../outputDisplay/PseudoCaloMET_Nm1_%srescaled_bis.pdf", OTHERrescale ? "OTHER" : ""));

    return;
}

void nHSCP(const char *qcd = "../output/QCD2024_V16/QCD2024_mu_V16p1_weighted.root",
           const char *ttbar = "../output/TTbar2024_V15/TTbar2024_V15p6_weighted.root",
           const char *wjets = "../output/Wjets2024_V14/WjetMuNu2024_V14p8_weighted.root",
           const char *data = "../output/JetMET2024_V12/JetMET2024_V12p24.root") {

    gErrorIgnoreLevel = kError;

    // backgrounds in stacking order; an empty path leaves the sample out
    struct Bkg { const char *path; const char *key; const char *label; int color; TFile *f; };
    std::vector<Bkg> bkgs;
    for (Bkg b : std::vector<Bkg>{{qcd, "qcd", "QCD muEnriched", kGreen-4, nullptr},
                                  {ttbar, "ttbar", "TTbar", kRed, nullptr},
                                  {wjets, "wjets", "W+jets", kBlue-7, nullptr}}) {
        b.f = OpenOptional(b.path);
        if (b.f) bkgs.push_back(b);
        else std::cout << "nHSCP: " << b.key << " left out" << std::endl;
    }
    TFile *ifileMET = new TFile(data, "READ");

    TH1F *nHSCP_MET = (TH1F*)ifileMET->Get("nHSCP");
    TH1F *Sel_nHSCP_MET = (TH1F*)ifileMET->Get("METanalysis_Eta2p4_nHSCP");

    TLegend *legend = new TLegend(0.7, 0.5, 0.85, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

    // data on a stack of the backgrounds, for one histogram (the same name in every file)
    auto drawStack = [&](TH1F *hData, const char *tag) {
        if (!hData) Abort(Form("%s: a histogram is missing from the data file %s", tag, data));
        TCanvas *c = new TCanvas(Form("c_%s", tag), Form("c_%s", tag), 800, 600);
        THStack *hs = new THStack(Form("hs_%s", tag), "");
        bool fillLegend = legend->GetListOfPrimitives()->GetSize() == 0;
        for (auto &b : bkgs) {
            TH1F *h = (TH1F*)b.f->Get(hData->GetName());
            if (!h) Abort(Form("nHSCP: no %s in %s", hData->GetName(), b.path));
            h->SetLineColor(b.color);
            h->SetFillColorAlpha(b.color, 0.5);
            hs->Add(h);
            if (fillLegend) legend->AddEntry(h, Label(b.key, b.label), "f");
        }
        hData->SetLineColor(kBlack);
        hData->SetMarkerColor(kBlack);
        hData->SetMarkerStyle(20);
        if (fillLegend) legend->AddEntry(hData, Label("data", "JetMET data"), "lep");
        TH1 *frame = hData;
        if (!bkgs.empty()) {
            hs->Draw("HIST");
            frame = hs->GetHistogram();
            hs->SetMinimum(0.8);
        }
        else hData->SetMinimum(0.8);
        frame->GetXaxis()->SetTitle("N_{HSCP}");
        frame->GetYaxis()->SetTitle("Number of events");
        hData->Draw(bkgs.empty() ? "E1" : "E1 same");
        legend->Draw();
        c->SetLogy();
        return c;
    };

    drawStack(nHSCP_MET, "nHSCP")->SaveAs(OutPath("nHSCP_nosel.pdf"));
    drawStack(Sel_nHSCP_MET, "Sel_nHSCP")->SaveAs(OutPath("nHSCP_sel.pdf"));

    return;
}

void nPV(const char *qcd = "../output/QCD2024_V16/QCD2024_mu_V16p1_weighted.root",
         const char *ttbar = "../output/TTbar2024_V15/TTbar2024_V15p6_weighted.root",
         const char *wjets = "../output/Wjets2024_V14/WjetMuNu2024_V14p8_weighted.root",
         const char *data = "../output/JetMET2024_V12/JetMET2024_V12p24.root") {

    gErrorIgnoreLevel = kError;

    // backgrounds in stacking order; an empty path leaves the sample out
    struct Bkg { const char *path; const char *key; const char *label; int color; TFile *f; };
    std::vector<Bkg> bkgs;
    for (Bkg b : std::vector<Bkg>{{qcd, "qcd", "QCD muEnriched", kGreen-4, nullptr},
                                  {ttbar, "ttbar", "TTbar", kRed, nullptr},
                                  {wjets, "wjets", "W+jets", kBlue-7, nullptr}}) {
        b.f = OpenOptional(b.path);
        if (b.f) bkgs.push_back(b);
        else std::cout << "nPV: " << b.key << " left out" << std::endl;
    }
    TFile *ifileMET = new TFile(data, "READ");

    TH1F *nPV_MET = (TH1F*)ifileMET->Get("nPVgood");
    TH1F *Sel_nPV_MET = (TH1F*)ifileMET->Get("METanalysis_Eta2p4_nPVgood");

    TLegend *legend = new TLegend(0.7, 0.5, 0.85, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

    // data on a stack of the backgrounds, for one histogram (the same name in every file)
    auto drawStack = [&](TH1F *hData, const char *tag) {
        if (!hData) Abort(Form("%s: a histogram is missing from the data file %s", tag, data));
        TCanvas *c = new TCanvas(Form("c_%s", tag), Form("c_%s", tag), 800, 600);
        THStack *hs = new THStack(Form("hs_%s", tag), "");
        bool fillLegend = legend->GetListOfPrimitives()->GetSize() == 0;
        for (auto &b : bkgs) {
            TH1F *h = (TH1F*)b.f->Get(hData->GetName());
            if (!h) Abort(Form("nPV: no %s in %s", hData->GetName(), b.path));
            h->SetLineColor(b.color);
            h->SetFillColorAlpha(b.color, 0.5);
            hs->Add(h);
            if (fillLegend) legend->AddEntry(h, Label(b.key, b.label), "f");
        }
        hData->SetLineColor(kBlack);
        hData->SetMarkerColor(kBlack);
        hData->SetMarkerStyle(20);
        if (fillLegend) legend->AddEntry(hData, Label("data", "JetMET data"), "lep");
        TH1 *frame = hData;
        if (!bkgs.empty()) {
            hs->Draw("HIST");
            frame = hs->GetHistogram();
            hs->SetMinimum(0.8);
        }
        else hData->SetMinimum(0.8);
        frame->GetXaxis()->SetTitle("N_{PV}");
        frame->GetYaxis()->SetTitle("Number of events");
        hData->Draw(bkgs.empty() ? "E1" : "E1 same");
        legend->Draw();
        c->SetLogy();
        return c;
    };

    drawStack(nPV_MET, "nPV")->SaveAs(OutPath("nPV_nosel.pdf"));
    drawStack(Sel_nPV_MET, "Sel_nPV")->SaveAs(OutPath("nPV_sel.pdf"));

    return;
}

void PostTriggerPseudoMET() {

    gErrorIgnoreLevel = kError;

    TFile *ifileWjetMuNu = new TFile("../output/Wjets2024_V14/WjetMuNu2024_V14p9_weighted.root", "READ");
    TFile *ifileMET = new TFile("../output/JetMET2024_V12/JetMET2024_V12p25.root", "READ");

    TH1F *PostTrigger_PseudoCaloMET_SF_rescaled_Wmunu = (TH1F*)ifileWjetMuNu->Get("PostTrigger_PseudoCaloMET_SF_rescaled");
    TH1F *PostTrigger_PseudoCaloMET_SF_NOTrescaled_Wmunu = (TH1F*)ifileWjetMuNu->Get("PostTrigger_PseudoCaloMET_SF_NOTrescaled");
    TH1F *PostTrigger_PseudoCaloMET_MET = (TH1F*)ifileMET->Get("PostTrigger_PseudoCaloMET_SF_NOTrescaled");

    PostTrigger_PseudoCaloMET_SF_rescaled_Wmunu->Scale(1./PostTrigger_PseudoCaloMET_SF_rescaled_Wmunu->Integral());
    PostTrigger_PseudoCaloMET_SF_NOTrescaled_Wmunu->Scale(1./PostTrigger_PseudoCaloMET_SF_NOTrescaled_Wmunu->Integral());
    PostTrigger_PseudoCaloMET_MET->Scale(1./PostTrigger_PseudoCaloMET_MET->Integral());

    PostTrigger_PseudoCaloMET_SF_rescaled_Wmunu->SetLineColor(kBlue-7);
    PostTrigger_PseudoCaloMET_SF_rescaled_Wmunu->SetFillColorAlpha(kBlue-7, 0.5);
    PostTrigger_PseudoCaloMET_SF_NOTrescaled_Wmunu->SetLineColor(kCyan);
    PostTrigger_PseudoCaloMET_SF_NOTrescaled_Wmunu->SetFillColorAlpha(kCyan, 0.5);
    PostTrigger_PseudoCaloMET_MET->SetLineColor(kBlack);
    PostTrigger_PseudoCaloMET_MET->SetMarkerColor(kBlack);
    PostTrigger_PseudoCaloMET_MET->SetMarkerStyle(20);

    TLegend *legend = new TLegend(0.7, 0.5, 0.85, 0.9);
    legend->AddEntry(PostTrigger_PseudoCaloMET_SF_rescaled_Wmunu, "WjetMuNu SFwRescaling", "f");
    legend->AddEntry(PostTrigger_PseudoCaloMET_SF_NOTrescaled_Wmunu, "WjetMuNu SFwo/Rescaling", "f");
    legend->AddEntry(PostTrigger_PseudoCaloMET_MET, "MET2024", "lep");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

    TCanvas *c_PseudoCaloMET = new TCanvas("c_PseudoCaloMET", "c_PseudoCaloMET", 800, 600);
    c_PseudoCaloMET->cd();
    PostTrigger_PseudoCaloMET_SF_rescaled_Wmunu->GetXaxis()->SetTitle("PseudoMET (~ CaloMET) [GeV]");
    PostTrigger_PseudoCaloMET_SF_rescaled_Wmunu->GetYaxis()->SetTitle("Number of events");
    PostTrigger_PseudoCaloMET_SF_rescaled_Wmunu->GetXaxis()->SetRangeUser(20, 1000);
    PostTrigger_PseudoCaloMET_SF_rescaled_Wmunu->Draw("HIST");
    PostTrigger_PseudoCaloMET_SF_NOTrescaled_Wmunu->Draw("HIST same");
    PostTrigger_PseudoCaloMET_MET->Draw("E1 same");
    legend->Draw();
    gPad->SetLogy();
    gStyle->SetOptStat(0);


    TCanvas *cRatio_Nm1 = DrawWithRatio(PostTrigger_PseudoCaloMET_MET, 
                                        PostTrigger_PseudoCaloMET_SF_rescaled_Wmunu,
                                        PostTrigger_PseudoCaloMET_SF_NOTrescaled_Wmunu,
                                        c_PseudoCaloMET,
                                        "PseudoMET w and wo rescaling (SF applied)",
                                        "data/MC",
                                        "PseudoMET (~ CaloMET) [GeV]",
                                        20,
                                        1000);

    cRatio_Nm1->SaveAs("../outputDisplay/PostTriggerLOG_PseudoCaloMET_SF_rescaled_vs_NOTrescaled.pdf");

    TCanvas *c_PseudoCaloMET_b = new TCanvas("c_PseudoCaloMET_b", "c_PseudoCaloMET_b", 800, 600);
    c_PseudoCaloMET_b->cd();
    PostTrigger_PseudoCaloMET_SF_rescaled_Wmunu->GetXaxis()->SetTitle("PseudoMET (~ CaloMET) [GeV]");
    PostTrigger_PseudoCaloMET_SF_rescaled_Wmunu->GetYaxis()->SetTitle("Number of events");
    PostTrigger_PseudoCaloMET_SF_rescaled_Wmunu->GetXaxis()->SetRangeUser(20, 1000);
    PostTrigger_PseudoCaloMET_SF_rescaled_Wmunu->SetMaximum(0.08);
    PostTrigger_PseudoCaloMET_SF_rescaled_Wmunu->Draw("HIST");
    PostTrigger_PseudoCaloMET_SF_NOTrescaled_Wmunu->Draw("HIST same");
    PostTrigger_PseudoCaloMET_MET->Draw("E1 same");
    legend->Draw();
    gStyle->SetOptStat(0);

    TCanvas *cRatio_Nm1_b = DrawWithRatio(PostTrigger_PseudoCaloMET_MET, 
                                        PostTrigger_PseudoCaloMET_SF_rescaled_Wmunu,
                                        PostTrigger_PseudoCaloMET_SF_NOTrescaled_Wmunu,
                                        c_PseudoCaloMET_b,
                                        "PseudoMET w and wo rescaling (SF applied)",
                                        "data/MC",
                                        "PseudoMET (~ CaloMET) [GeV]",
                                        20,
                                        1000);

    cRatio_Nm1_b->SaveAs("../outputDisplay/PostTrigger_PseudoCaloMET_SF_rescaled_vs_NOTrescaled.pdf");

    return;
}

void MyClusters() {

    // build a histogram of 512 bins with x axis from 0 to 512, where two strip clusters are artificially built:
    TH1F *h_clusters = new TH1F("", "", 512, -0.5, 511.5);
    h_clusters->SetBinContent(103, 6);
    h_clusters->SetBinContent(104, 18);
    h_clusters->SetBinContent(105, 177);
    h_clusters->SetBinContent(106, 22);
    h_clusters->SetBinContent(107, 7);

    h_clusters->SetBinContent(436, 4);
    h_clusters->SetBinContent(437, 129);
    h_clusters->SetBinContent(438, 83);
    h_clusters->SetBinContent(439, 15);
    h_clusters->SetBinContent(440, 2);

    h_clusters->GetXaxis()->SetTitle("Strip index");
    h_clusters->GetYaxis()->SetTitle("ADC counts");
    h_clusters->GetYaxis()->SetTitleSize(0.06);
    h_clusters->GetXaxis()->SetTitleSize(0.06);
    h_clusters->GetXaxis()->SetTitleOffset(0.9);
    h_clusters->GetYaxis()->SetTitleOffset(1.3);
    h_clusters->GetXaxis()->SetLabelSize(0.05);
    h_clusters->GetYaxis()->SetLabelSize(0.05);
    h_clusters->SetLineColor(kBlue-7);
    h_clusters->SetLineWidth(2);
    h_clusters->SetFillColorAlpha(kBlue-7, 0.5);
    gStyle->SetOptStat(0);

    TLatex *latex1 = new TLatex(0.155, 0.91, "#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);

    // clone of the histo
    TH1F *h_clusters_clone = (TH1F*)h_clusters->Clone("h_clusters_clone");

    h_clusters_clone->GetYaxis()->SetTitleSize(0.06);
    h_clusters_clone->GetXaxis()->SetTitleSize(0.06);
    h_clusters_clone->GetXaxis()->SetTitleOffset(0.9);
    h_clusters_clone->GetYaxis()->SetTitleOffset(1.3);
    h_clusters_clone->GetXaxis()->SetLabelSize(0.05);
    h_clusters_clone->GetYaxis()->SetLabelSize(0.05);

    TCanvas *c_clusters1 = new TCanvas("c_clusters1", "c_clusters1", 800, 800);
    c_clusters1->cd();
    c_clusters1->SetLeftMargin(0.16); c_clusters1->SetBottomMargin(0.16);
    h_clusters->GetXaxis()->SetRangeUser(101, 107);
    h_clusters->Draw("HIST");
    latex1->Draw("same");
    c_clusters1->SaveAs("../outputDisplay/Cluster_centered.pdf");
    latex1->SetTitle("#it{Private work (CMS simulation)}");
    c_clusters1->Modified();
    c_clusters1->Update();
    c_clusters1->SaveAs("../outputDisplay/Cluster_centered_bis.pdf");

    TCanvas *c_clusters2 = new TCanvas("c_clusters2", "c_clusters2", 800, 800);
    c_clusters2->cd();
    c_clusters2->SetLeftMargin(0.16); c_clusters2->SetBottomMargin(0.16);
    h_clusters_clone->GetXaxis()->SetRangeUser(434, 440);
    h_clusters_clone->Draw("HIST");
    latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    latex1->Draw("same");
    c_clusters2->SaveAs("../outputDisplay/Cluster_tilted.pdf");
    latex1->SetTitle("#it{Private work (CMS simulation)}");
    c_clusters2->Modified();
    c_clusters2->Update();
    c_clusters2->SaveAs("../outputDisplay/Cluster_tilted_bis.pdf");


    TH1F *hsat = new TH1F("", "", 512, -0.5, 511.5);
    hsat->SetBinContent(211, 16);
    hsat->SetBinContent(212, 38);
    hsat->SetBinContent(213, 391);
    hsat->SetBinContent(214, 43);
    hsat->SetBinContent(215, 18);
    hsat->SetBinContent(216, 3);

    hsat->GetYaxis()->SetTitleSize(0.06);
    hsat->GetXaxis()->SetTitleSize(0.06);
    hsat->GetXaxis()->SetTitleOffset(0.9);
    hsat->GetYaxis()->SetTitleOffset(1.3);
    hsat->GetXaxis()->SetLabelSize(0.05);
    hsat->GetYaxis()->SetLabelSize(0.05);


    TH1F *hsatreco = new TH1F("", "", 512, -0.5, 511.5);
    hsatreco->SetBinContent(211, 14);
    hsatreco->SetBinContent(212, 40);
    hsatreco->SetBinContent(213, 254);
    hsatreco->SetBinContent(214, 41);
    hsatreco->SetBinContent(215, 18);

    hsat->GetXaxis()->SetTitle("Strip index");
    hsat->GetYaxis()->SetTitle("ADC counts");
    hsat->SetLineColor(kBlue-7);
    hsat->SetLineWidth(2);
    hsat->SetFillColorAlpha(kBlue-7, 0.5);
    hsatreco->SetLineColor(kGreen-4);
    hsatreco->SetLineWidth(3);
    hsatreco->SetMarkerColor(kGreen-4);
    hsatreco->SetMarkerStyle(22);
    hsatreco->SetMarkerSize(2);
    gStyle->SetOptStat(0);

    TCanvas *c_sat = new TCanvas("c_sat", "c_sat", 800, 800);
    c_sat->cd();
    c_sat->SetLeftMargin(0.16); c_sat->SetBottomMargin(0.16);
    hsat->GetXaxis()->SetRangeUser(209, 216);
    hsat->Draw("HIST");
    hsatreco->Draw("hist same");
    hsatreco->Draw("P same");
    TLegend *legend = new TLegend(0.55, 0.75, 0.95, 0.9);
    legend->AddEntry(hsat, "Simulated", "f");
    legend->AddEntry(hsatreco, "Reconstructed", "lp");
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->Draw();
    latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    latex1->Draw("same");
    c_sat->SaveAs("../outputDisplay/Cluster_saturation.pdf");
    latex1->SetTitle("#it{Private work (CMS simulation)}");
    c_sat->Modified();
    c_sat->Update();
    c_sat->SaveAs("../outputDisplay/Cluster_saturation_bis.pdf");


    TH1F *h_LR1 = new TH1F("", "", 512, -0.5, 511.5);
    h_LR1->SetBinContent(213, 10);
    h_LR1->SetBinContent(214, 187);
    h_LR1->SetBinContent(215, 351);
    h_LR1->SetBinContent(216, 221);
    h_LR1->SetBinContent(217, 11);
    h_LR1->SetBinContent(218, 2);

    TH1F *h_LR1reco = new TH1F("", "", 512, -0.5, 511.5);
    h_LR1reco->SetBinContent(213, 13);
    h_LR1reco->SetBinContent(214, 160);
    h_LR1reco->SetBinContent(215, 254);
    h_LR1reco->SetBinContent(216, 172);
    h_LR1reco->SetBinContent(217, 12);

    h_LR1->SetLineColor(kBlue-7);
    h_LR1->SetLineWidth(2);
    h_LR1->SetFillColorAlpha(kBlue-7, 0.5);
    h_LR1reco->SetLineColor(kGreen-4);
    h_LR1reco->SetLineWidth(3);
    h_LR1reco->SetMarkerColor(kGreen-4);
    h_LR1reco->SetMarkerStyle(22);
    h_LR1reco->SetMarkerSize(2);

    TH1F *h_LR2 = new TH1F("", "", 512, -0.5, 511.5);
    h_LR2->SetBinContent(85, 2);
    h_LR2->SetBinContent(86, 14);
    h_LR2->SetBinContent(87, 159);
    h_LR2->SetBinContent(88, 271);
    h_LR2->SetBinContent(89, 231);
    h_LR2->SetBinContent(90, 15);
    h_LR2->SetBinContent(91, 4);

    TH1F *h_LR2reco = new TH1F("", "", 512, -0.5, 511.5);
    h_LR2reco->SetBinContent(87, 148);
    h_LR2reco->SetBinContent(88, 254);
    h_LR2reco->SetBinContent(89, 225);
    h_LR2reco->SetBinContent(90, 17);

    h_LR2->SetLineColor(kBlue-7);
    h_LR2->SetLineWidth(2);
    h_LR2->SetFillColorAlpha(kBlue-7, 0.5);
    h_LR2reco->SetLineColor(kGreen-4);
    h_LR2reco->SetLineWidth(3);
    h_LR2reco->SetMarkerColor(kGreen-4);
    h_LR2reco->SetMarkerStyle(22);
    h_LR2reco->SetMarkerSize(2);

    TH1F *h_LR3 = new TH1F("", "", 512, -0.5, 511.5);
    h_LR3->SetBinContent(478, 2);
    h_LR3->SetBinContent(479, 14);
    h_LR3->SetBinContent(480, 174);
    h_LR3->SetBinContent(481, 241);
    h_LR3->SetBinContent(482, 289);
    h_LR3->SetBinContent(483, 222);
    h_LR3->SetBinContent(484, 15);
    h_LR3->SetBinContent(485, 3);

    TH1F *h_LR3reco = new TH1F("", "", 512, -0.5, 511.5);
    h_LR3reco->SetBinContent(480, 168);
    h_LR3reco->SetBinContent(481, 219);
    h_LR3reco->SetBinContent(482, 254);
    h_LR3reco->SetBinContent(483, 198);
    h_LR3reco->SetBinContent(484, 35);

    h_LR3->SetLineColor(kBlue-7);
    h_LR3->SetLineWidth(2);
    h_LR3->SetFillColorAlpha(kBlue-7, 0.5);
    h_LR3reco->SetLineColor(kGreen-4);
    h_LR3reco->SetLineWidth(3);
    h_LR3reco->SetMarkerColor(kGreen-4);
    h_LR3reco->SetMarkerStyle(22);
    h_LR3reco->SetMarkerSize(2);

    h_LR1->GetYaxis()->SetTitleSize(0.06);
    h_LR1->GetXaxis()->SetTitleSize(0.06);
    h_LR1->GetXaxis()->SetTitleOffset(0.9);
    h_LR1->GetYaxis()->SetTitleOffset(1.3);
    h_LR1->GetXaxis()->SetLabelSize(0.05);
    h_LR1->GetYaxis()->SetLabelSize(0.05);

    h_LR2->GetYaxis()->SetTitleSize(0.06);
    h_LR2->GetXaxis()->SetTitleSize(0.06);
    h_LR2->GetXaxis()->SetTitleOffset(0.9);
    h_LR2->GetYaxis()->SetTitleOffset(1.3);
    h_LR2->GetXaxis()->SetLabelSize(0.04);
    h_LR2->GetYaxis()->SetLabelSize(0.04);

    h_LR3->GetYaxis()->SetTitleSize(0.06);
    h_LR3->GetXaxis()->SetTitleSize(0.06);
    h_LR3->GetXaxis()->SetTitleOffset(0.9);
    h_LR3->GetYaxis()->SetTitleOffset(1.3);
    h_LR3->GetXaxis()->SetLabelSize(0.04);
    h_LR3->GetYaxis()->SetLabelSize(0.04);

    

    TCanvas *c_1 = new TCanvas("c_1", "c_1", 800, 800);
    c_1->cd();
    c_1->SetLeftMargin(0.16); c_1->SetBottomMargin(0.16);
    h_LR1->GetXaxis()->SetRangeUser(211, 218);
    h_LR1->Draw("HIST");
    h_LR1->GetXaxis()->SetTitle("Strip index");
    h_LR1->GetYaxis()->SetTitle("ADC counts");
    h_LR1reco->Draw("hist same");
    h_LR1reco->Draw("P same");
    TLegend *legend2 = new TLegend(0.18, 0.74, 0.4, 0.9);
    legend2->AddEntry(hsat, "Simulated", "f");
    legend2->AddEntry(hsatreco, "Reconstructed", "lp");
    legend2->SetTextSize(0.035);
    legend2->SetBorderSize(0);
    legend2->SetFillStyle(0);
    legend2->Draw();
    latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    latex1->Draw("same");
    c_1->SaveAs("../outputDisplay/LR_overcorrected_1.pdf");
    latex1->SetTitle("#it{Private work (CMS simulation)}");
    c_1->Modified();
    c_1->Update();
    c_1->SaveAs("../outputDisplay/LR_overcorrected_1_bis.pdf");

    TCanvas *c_2 = new TCanvas("c_2", "c_2", 800, 800);
    c_2->cd();
    c_2->SetLeftMargin(0.16); c_2->SetBottomMargin(0.16);
    h_LR2->GetXaxis()->SetRangeUser(83, 91);
    h_LR2->Draw("HIST");
    h_LR2->GetXaxis()->SetTitle("Strip index");
    h_LR2->GetYaxis()->SetTitle("ADC counts");
    h_LR2reco->Draw("hist same");
    h_LR2reco->Draw("P same");
    legend2->Draw();
    latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    latex1->Draw("same");
    c_2->SaveAs("../outputDisplay/LR_overcorrected_2.pdf");
    latex1->SetTitle("#it{Private work (CMS simulation)}");
    c_2->Modified();
    c_2->Update();
    c_2->SaveAs("../outputDisplay/LR_overcorrected_2_bis.pdf");

    TCanvas *c_3 = new TCanvas("c_3", "c_3", 800, 800);
    c_3->cd();
    c_3->SetLeftMargin(0.16); c_3->SetBottomMargin(0.16);
    h_LR3->GetXaxis()->SetRangeUser(476, 485);
    h_LR3->Draw("HIST");
    h_LR3->GetXaxis()->SetTitle("Strip index");
    h_LR3->GetYaxis()->SetTitle("ADC counts");
    h_LR3reco->Draw("hist same");
    h_LR3reco->Draw("P same");
    legend2->Draw();
    latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    latex1->Draw("same");
    c_3->SaveAs("../outputDisplay/LR_overcorrected_3.pdf");
    latex1->SetTitle("#it{Private work (CMS simulation)}");
    c_3->Modified();
    c_3->Update();
    c_3->SaveAs("../outputDisplay/LR_overcorrected_3_bis.pdf");


    return;
}

void Run2_vs_Run3_gluino__TriggerEff() {

    TFile *ifileRun2 = new TFile("../output/Gluino_V19/Gluino_Run2_MET_madgraph_2000_V13p2.root", "READ");
    TFile *ifileRun3 = new TFile("../output/Gluino_V19/Gluino_Run3_MET_madgraph_2000_V19p1.root", "READ");

    TH1F *PseudoCaloMET_Run2 = (TH1F*)ifileRun2->Get("METanalysis_Eta2p4_EffTrg_PseudoCaloMET");
    TH1F *if___orMETtrg___PseudoCaloMET_Run2 = (TH1F*)ifileRun2->Get("METanalysis_Eta2p4_EffTrg_if___orMETtrg___PseudoCaloMET");
    TH1F *PseudoCaloMET_Run3 = (TH1F*)ifileRun3->Get("METanalysis_Eta2p4_EffTrg_PseudoCaloMET");
    TH1F *if___orMETtrg___PseudoCaloMET_Run3 = (TH1F*)ifileRun3->Get("METanalysis_Eta2p4_EffTrg_if___orMETtrg___PseudoCaloMET");
    
    TH1F *eff_orMETtrg_PseudoCaloMET_Run2 = (TH1F*)if___orMETtrg___PseudoCaloMET_Run2->Clone("eff_orMETtrg_PseudoCaloMET_Run2");
    TH1F *eff_orMETtrg_PseudoCaloMET_Run3 = (TH1F*)if___orMETtrg___PseudoCaloMET_Run3->Clone("eff_orMETtrg_PseudoCaloMET_Run3");

    eff_orMETtrg_PseudoCaloMET_Run2->Divide(PseudoCaloMET_Run2);
    eff_orMETtrg_PseudoCaloMET_Run3->Divide(PseudoCaloMET_Run3);

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    gStyle->SetOptStat(0);
    eff_orMETtrg_PseudoCaloMET_Run2->GetXaxis()->SetTitle("PseudoMET (~ CaloMET) [GeV]");
    eff_orMETtrg_PseudoCaloMET_Run2->GetYaxis()->SetTitle("eff. OR all MET triggers");
    eff_orMETtrg_PseudoCaloMET_Run2->GetYaxis()->SetRangeUser(0, 1);
    eff_orMETtrg_PseudoCaloMET_Run2->GetXaxis()->SetRangeUser(0, 1200);


    eff_orMETtrg_PseudoCaloMET_Run2->SetLineColor(kRed);
    eff_orMETtrg_PseudoCaloMET_Run2->SetMarkerColor(kRed);
    eff_orMETtrg_PseudoCaloMET_Run2->SetMarkerStyle(20);
    eff_orMETtrg_PseudoCaloMET_Run2->Draw("E1");
    eff_orMETtrg_PseudoCaloMET_Run3->SetLineColor(kBlue-7);
    eff_orMETtrg_PseudoCaloMET_Run3->SetMarkerColor(kBlue-7);
    eff_orMETtrg_PseudoCaloMET_Run3->SetMarkerStyle(22);
    eff_orMETtrg_PseudoCaloMET_Run3->Draw("E1 same");

    TLegend *legend = new TLegend(0.6, 0.2, 0.85, 0.5);
    legend->AddEntry(eff_orMETtrg_PseudoCaloMET_Run2, "Run2", "lep");
    legend->AddEntry(eff_orMETtrg_PseudoCaloMET_Run3, "Run3", "lep");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->Draw();

    c1->SaveAs("../outputDisplay/Run2_vs_Run3_gluino__TriggerEff.pdf");

    return;
}

void ComparePseudoMET(bool isRescaled = true) {

    TFile *ifileWmunu = new TFile("../output/Wjets2024_V14/WjetMuNu2024_V14p9_weighted.root", "READ");
    TFile *ifileTTbar = new TFile("../output/TTbar2024_V15/TTbar2024_V15p8.root", "READ");
    TFile *ifileGluino = new TFile("../output/Gluino_V19/Gluino_Run3_MET_madgraph_2000_V19p1.root", "READ");

    TH1F *PseudoMET_Wmunu = (TH1F*)ifileWmunu->Get(Form("Nosel_PseudoCaloMET%s", isRescaled ? "_weighted" : ""));
    TH1F *PseudoMET_TTbar = (TH1F*)ifileTTbar->Get(Form("Nosel_PseudoCaloMET%s", isRescaled ? "_weighted" : ""));
    TH1F *PseudoMET_Gluino = (TH1F*)ifileGluino->Get(Form("Nosel_PseudoCaloMET%s", isRescaled ? "_rescaled" : ""));

    PseudoMET_Wmunu->Scale(1./PseudoMET_Wmunu->Integral());
    PseudoMET_TTbar->Scale(1./PseudoMET_TTbar->Integral());
    PseudoMET_Gluino->Scale(1./PseudoMET_Gluino->Integral());

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    gStyle->SetOptStat(0);
    PseudoMET_Wmunu->GetXaxis()->SetTitle("PseudoMET (~ CaloMET) [GeV]");
    PseudoMET_Wmunu->GetYaxis()->SetTitle("Events (normalised)");

    PseudoMET_Wmunu->SetLineColor(kBlue-7);
    PseudoMET_Wmunu->SetFillColorAlpha(kBlue-7, 0.5);
    PseudoMET_Wmunu->Draw("HIST");
    PseudoMET_TTbar->SetLineColor(kRed);
    PseudoMET_TTbar->SetFillColorAlpha(kRed, 0.5);
    PseudoMET_TTbar->Draw("HIST same");
    PseudoMET_Gluino->SetLineColor(kBlack);
    PseudoMET_Gluino->Draw("HIST same");
    TLegend *legend = new TLegend(0.6, 0.7, 0.85, 0.9);
    legend->AddEntry(PseudoMET_Wmunu, "W+jets", "f");
    legend->AddEntry(PseudoMET_TTbar, "TTbar", "f");
    legend->AddEntry(PseudoMET_Gluino, "Gluino", "f");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->Draw();
    c1->SetLogy();
    c1->SaveAs(Form("../outputDisplay/Compare_PseudoMET%s.pdf", isRescaled ? "_rescaled" : "_NOTrescaled"));

    return;
}

void YieldAfterSF(bool isRescaled = true, bool isMuSelection = false, bool isMETSelection = false, bool isTTbarSelection = false,
                  const char *wjets = "../output/Wjets2024_V14/WjetMuNu2024_V14p8_weighted.root",
                  const char *ttbar = "../output/TTbar2024_V15/TTbar2024_V15p8_weighted.root",
                  const char *dataMu = "../output/Mu2024_V18/Mu2024_V18.root",
                  const char *dataMET = "../output/JetMET2024_V12/JetMET2024_V12p24.root",
                  const char *dataMuonEG = "../output/MuonEG_V17/MuonEG2024_V17p4.root",
                  const char *gluino = "../output/Gluino_V19/Gluino_Run3_MET_madgraph_2000_V19p1.root") {

    // only the files of the chosen selection are needed; the others may be empty paths
    TFile *ifileWmunu = OpenOptional(wjets);
    TFile *ifileTTbar = OpenOptional(ttbar);
    TFile *ifileDataMu = OpenOptional(dataMu);
    TFile *ifileDataMET = OpenOptional(dataMET);
    TFile *ifileDataMuonEG = OpenOptional(dataMuonEG);
    TFile *ifileGluino = OpenOptional(gluino);
    auto need = [](TFile *f, const char *what) {
        if (!f) Abort(Form("YieldAfterSF: this selection needs %s", what));
        return f;
    };

    TH1F *PostTrigger_PseudoCaloMET_SF_rescaled___MC, *PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA, *PostTrigger_PseudoCaloMET_SF_rescaled___Gluino;
    if (isMuSelection) {
        PostTrigger_PseudoCaloMET_SF_rescaled___MC = (TH1F*)need(ifileWmunu, "wjets")->Get(Form("CalibPseudoMET_MuWay%s_if___orMETtrg___PseudoCaloMET", isRescaled ? "_isRescaled" : ""));
        PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA = (TH1F*)need(ifileDataMu, "dataMu")->Get("CalibPseudoMET_MuWay_if___orMETtrg___PseudoCaloMET"); // take the not rescaled one for data.
        PostTrigger_PseudoCaloMET_SF_rescaled___Gluino = (TH1F*)need(ifileGluino, "gluino")->Get(Form("METanalysis_Eta2p4_PseudoCaloMET%s", isRescaled ? "_rescaled" : ""));
    }
    else if (isMETSelection) {
        PostTrigger_PseudoCaloMET_SF_rescaled___MC = (TH1F*)need(ifileWmunu, "wjets")->Get(Form("Nm1_event_CaloMET%s", isRescaled ? "_weighted" : ""));
        PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA = (TH1F*)need(ifileDataMET, "dataMET")->Get("Nm1_event_CaloMET"); // take the not rescaled one for data.
        PostTrigger_PseudoCaloMET_SF_rescaled___Gluino = (TH1F*)need(ifileGluino, "gluino")->Get(Form("Nm1_event_CaloMET%s", isRescaled ? "_rescaled" : ""));
    }
    else if (isTTbarSelection) {
        PostTrigger_PseudoCaloMET_SF_rescaled___MC = (TH1F*)need(ifileTTbar, "ttbar")->Get(Form("CalibPseudoMET%s_if___orMETtrg___PseudoCaloMET", isRescaled ? "_isRescaled" : ""));
        PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA = (TH1F*)need(ifileDataMuonEG, "dataMuonEG")->Get("CalibPseudoMET_if___orMETtrg___PseudoCaloMET"); // take the not rescaled one for data.
        PostTrigger_PseudoCaloMET_SF_rescaled___Gluino = (TH1F*)need(ifileGluino, "gluino")->Get(Form("METanalysis_Eta2p4_PseudoCaloMET%s", isRescaled ? "_rescaled" : ""));
    }


    PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA->Scale(1./PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA->Integral());
    PostTrigger_PseudoCaloMET_SF_rescaled___MC->Scale(1./PostTrigger_PseudoCaloMET_SF_rescaled___MC->Integral());
    PostTrigger_PseudoCaloMET_SF_rescaled___Gluino->Scale(1./PostTrigger_PseudoCaloMET_SF_rescaled___Gluino->Integral());

    float yieldMC = 0;
    int binStart = PostTrigger_PseudoCaloMET_SF_rescaled___MC->FindBin(170.0);
    int binEnd = PostTrigger_PseudoCaloMET_SF_rescaled___MC->GetNbinsX() + 1; // +1 = overflow
    for (int i = binStart; i <= binEnd; i++) yieldMC += PostTrigger_PseudoCaloMET_SF_rescaled___MC->GetBinContent(i);
    
    float yieldData = 0;
    binStart = PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA->FindBin(170.0);
    binEnd = PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA->GetNbinsX() + 1; // +1 = overflow
    for (int i = binStart; i <= binEnd; i++) yieldData += PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA->GetBinContent(i);

    float yieldGluino = 0;
    binStart = PostTrigger_PseudoCaloMET_SF_rescaled___Gluino->FindBin(170.0);
    binEnd = PostTrigger_PseudoCaloMET_SF_rescaled___Gluino->GetNbinsX() + 1; // +1 = overflow
    for (int i = binStart; i <= binEnd; i++) yieldGluino += PostTrigger_PseudoCaloMET_SF_rescaled___Gluino->GetBinContent(i);

    cout << (isRescaled ? "With rescaling:" : "Without rescaling:") << endl;
    if (isMETSelection) cout << "HLT MET + HSCP selection:" << endl;
    else if (isMuSelection) cout << "HLT mu + HSCP selection:" << endl;
    else if (isTTbarSelection) cout << "HLT MuonEG + HSCP selection:" << endl;
    cout << "Yield in data (not rescaled) = " << yieldData << endl;
    cout << "Yield in MC (rescaled) = " << yieldMC << endl;
    if (isMETSelection) cout << "Yield in Gluino MC (rescaled) = " << yieldGluino << endl;
    cout << "ratio (Data/MC)= " << yieldData/yieldMC << endl;

    // clone PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA only for data from 170 up to N+1 and draw hit hashed:
    TH1F *PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA_clone = (TH1F*)PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA->Clone("PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA_clone");
    PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA_clone->Reset();
    for (int i = PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA->FindBin(170.0); i <= PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA->GetNbinsX() + 1; i++) 
        PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA_clone->SetBinContent(i, PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA->GetBinContent(i));
    TH1F *PostTrigger_PseudoCaloMET_SF_rescaled___MC_clone = (TH1F*)PostTrigger_PseudoCaloMET_SF_rescaled___MC->Clone("PostTrigger_PseudoCaloMET_SF_rescaled___MC_clone");
    PostTrigger_PseudoCaloMET_SF_rescaled___MC_clone->Reset();
    for (int i = PostTrigger_PseudoCaloMET_SF_rescaled___MC->FindBin(170.0); i <= PostTrigger_PseudoCaloMET_SF_rescaled___MC->GetNbinsX() + 1; i++) 
        PostTrigger_PseudoCaloMET_SF_rescaled___MC_clone->SetBinContent(i, PostTrigger_PseudoCaloMET_SF_rescaled___MC->GetBinContent(i));
    TH1F *PostTrigger_PseudoCaloMET_SF_rescaled___Gluino_clone = (TH1F*)PostTrigger_PseudoCaloMET_SF_rescaled___Gluino->Clone("PostTrigger_PseudoCaloMET_SF_rescaled___Gluino_clone");
    PostTrigger_PseudoCaloMET_SF_rescaled___Gluino_clone->Reset();
    for (int i = PostTrigger_PseudoCaloMET_SF_rescaled___Gluino->FindBin(170.0); i <= PostTrigger_PseudoCaloMET_SF_rescaled___Gluino->GetNbinsX() + 1; i++) 
        PostTrigger_PseudoCaloMET_SF_rescaled___Gluino_clone->SetBinContent(i, PostTrigger_PseudoCaloMET_SF_rescaled___Gluino->GetBinContent(i));

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    gStyle->SetOptStat(0);
    PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA->SetLineColor(kBlack);
    PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA->SetMarkerColor(kBlack);
    PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA->SetMarkerStyle(20);
    PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA_clone->SetLineColor(kBlack);
    PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA_clone->SetFillColorAlpha(kBlack, 0.5);
    PostTrigger_PseudoCaloMET_SF_rescaled___MC->SetLineColor(kBlue-7);
    PostTrigger_PseudoCaloMET_SF_rescaled___MC->SetFillColorAlpha(kBlue-7, 0.5);
    PostTrigger_PseudoCaloMET_SF_rescaled___MC_clone->SetLineColor(kBlue-7);
    PostTrigger_PseudoCaloMET_SF_rescaled___MC_clone->SetFillColorAlpha(kBlue-7, 0.5);
    PostTrigger_PseudoCaloMET_SF_rescaled___Gluino_clone->SetLineColor(kRed);
    PostTrigger_PseudoCaloMET_SF_rescaled___Gluino_clone->SetFillColorAlpha(kRed, 0.5);
    PostTrigger_PseudoCaloMET_SF_rescaled___Gluino->SetLineColor(kRed);
    PostTrigger_PseudoCaloMET_SF_rescaled___Gluino->SetMarkerColor(kRed);
    PostTrigger_PseudoCaloMET_SF_rescaled___Gluino->SetMarkerStyle(22);

    PostTrigger_PseudoCaloMET_SF_rescaled___MC->GetXaxis()->SetTitle("PseudoMET (~ CaloMET) [GeV]");
    PostTrigger_PseudoCaloMET_SF_rescaled___MC->GetYaxis()->SetTitle("Events (normalised)");

    PostTrigger_PseudoCaloMET_SF_rescaled___MC->Draw("HIST");
    //PostTrigger_PseudoCaloMET_SF_rescaled___MC_clone->Draw("HIST same");
    //PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA_clone->Draw("HIST same");
    PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA->Draw("E1 same");
    if (isMETSelection) {
        //PostTrigger_PseudoCaloMET_SF_rescaled___Gluino_clone->Draw("HIST same");
        PostTrigger_PseudoCaloMET_SF_rescaled___Gluino->Draw("E1 same");
    }
    TLegend *legend = new TLegend(0.6, 0.7, 0.85, 0.9);
    legend->AddEntry(PostTrigger_PseudoCaloMET_SF_NOTrescaled___DATA, "Data", "lep");
    legend->AddEntry(PostTrigger_PseudoCaloMET_SF_rescaled___MC, Form("MC (%sSF)",(isRescaled) ? "rescaled + " : ""), "f");
    if (isMETSelection) legend->AddEntry(PostTrigger_PseudoCaloMET_SF_rescaled___Gluino, Form("%s (%sSF)", Label("gluino", "Gluino2000"), (isRescaled) ? "rescaled + " : ""), "lep");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->Draw();
    c1->SetLogy();
    if (isMETSelection) c1->SaveAs(OutPath("YieldAfterSF_HLTmet_%s.pdf", isRescaled ? "_rescaled" : "_NOTrescaled"));
    else if (isMuSelection) c1->SaveAs(OutPath("YieldAfterSF_HLTmu_%s.pdf", isRescaled ? "_rescaled" : "_NOTrescaled"));
    else if (isTTbarSelection) c1->SaveAs(OutPath("YieldAfterSF_HLTMuonEG_%s.pdf", isRescaled ? "_rescaled" : "_NOTrescaled"));


    return;
}

void PlotSF(const std::string& fileWmunu = "../outputDisplay/SF_TriggerEff_Mu2024_WMuNu_PseudoMETrescaled.txt",
           const std::string& fileTTbar  = "../outputDisplay/SF_TriggerEff_MuonEG2024_TTbar_PseudoMETrescaled.txt") {

  auto parse = [&](const std::string& fname,
                    const std::string& name,
                    const std::string& title) -> TH1F* {
    std::ifstream f(fname);
    if (!f.is_open()) { std::cerr << "Cannot open: " << fname << "\n"; return nullptr; }

    std::vector<double> bins, vals, errs;
    double b, vd, v, vs;
    while (f >> b >> vd >> v >> vs) {
      bins.push_back(b); vals.push_back(v); errs.push_back((vs - vd) / 2.0);
    }

    int    n    = bins.size();
    double step = (n > 1) ? (bins[1] - bins[0]) : 25.0;
    std::vector<double> edges;
    for (int i = 0; i < n; i++) edges.push_back(bins[i] - step/2.0);
    edges.push_back(bins[n-1] + step/2.0);

    TH1F* h = new TH1F(name.c_str(), title.c_str(), n, edges.data());
    for (int i = 0; i < n; i++) { h->SetBinContent(i+1, vals[i]); h->SetBinError(i+1, errs[i]); }
    return h;
  };

  TH1F* hW = parse(fileWmunu, "hWmunu", "Scale Factors;PseudoMET (~ CaloMET) [GeV];SF");
  TH1F* hT = parse(fileTTbar,  "hTTbar", "Scale Factors;PseudoMET (~ CaloMET) [GeV];SF");
  if (!hW || !hT) return;

  hW->SetLineColor(kBlue+1);  hW->SetMarkerColor(kBlue+1);  hW->SetMarkerStyle(20);
  hT->SetLineColor(kRed+1);   hT->SetMarkerColor(kRed+1);   hT->SetMarkerStyle(21);

  TCanvas* c = new TCanvas("cSF", "Scale Factors", 800, 600);
  c->SetLeftMargin(0.12); c->SetBottomMargin(0.12);
  gStyle->SetOptStat(0);

  hW->GetYaxis()->SetRangeUser(0.65, 1.10);
  hW->Draw("E1"); hT->Draw("E1 SAME");

  TLine* line = new TLine(hW->GetXaxis()->GetXmin(), 1.0,
                            hW->GetXaxis()->GetXmax(), 1.0);
  line->SetLineStyle(2); line->SetLineColor(kGray+2); line->Draw();

  TLegend* leg = new TLegend(0.65, 0.15, 0.88, 0.32);
  leg->AddEntry(hW, "SF Muon/Wmunu", "lep");
  leg->AddEntry(hT, "SF MuonEG/TTbar", "lep");
  leg->SetBorderSize(0); leg->Draw();

  c->SaveAs("../outputDisplay/Plot_SF___MuonWmunu_vs_MuonEGttbar.pdf");
}


void TableGluino(bool isRescaled) {

    TFile *ifileGluino = new TFile("Gluino_Run3_MET_madgraph_2000_V19p1.root", "READ"); // ../output/Gluino_V19/

    // Départ
    TH1F *Nosel_PseudoCaloMET = (TH1F*)ifileGluino->Get(Form("Nosel_PseudoCaloMET%s", isRescaled ? "_rescaled" : ""));

    // Trigger
    TH1F *PostTrigger_PseudoCaloMET_rescaled = (TH1F*)ifileGluino->Get(Form("PostTrigger_PseudoCaloMET_%srescaled", isRescaled ? "" : "NOT"));

    // Trigger + SF
    TH1F *PostTrigger_PseudoCaloMET_SF_rescaled___Gluino = (TH1F*)ifileGluino->Get(Form("PostTrigger_PseudoCaloMET_SF_%srescaled", isRescaled ? "" : "NOT"));

    // PseudoMET > 170
    TH1F *PostTrigger_PseudoCaloMET_SF_rescaled___Gluino_clone = (TH1F*)PostTrigger_PseudoCaloMET_SF_rescaled___Gluino->Clone("PostTrigger_PseudoCaloMET_SF_rescaled___Gluino_clone");
    PostTrigger_PseudoCaloMET_SF_rescaled___Gluino_clone->Reset();
    for (int i = PostTrigger_PseudoCaloMET_SF_rescaled___Gluino->FindBin(170.0); i <= PostTrigger_PseudoCaloMET_SF_rescaled___Gluino->GetNbinsX() + 1; i++) 
        PostTrigger_PseudoCaloMET_SF_rescaled___Gluino_clone->SetBinContent(i, PostTrigger_PseudoCaloMET_SF_rescaled___Gluino->GetBinContent(i));

    // HSCP sel
    TH1F *METanalysis_Eta2p4_PseudoCaloMET = (TH1F*)ifileGluino->Get(Form("METanalysis_PseudoMET%sescaled_Eta2p4_PseudoCaloMET_%srescaled", isRescaled ? "r" : "notR", isRescaled ? "" : "NOT"));


    cout << "# Gluino2000        " << (isRescaled ? "Rescaled" : "NOT Rescaled") << endl;
    cout << "Initial yield:       " << Nosel_PseudoCaloMET->Integral() << endl;
    cout << "After trigger:       " << PostTrigger_PseudoCaloMET_rescaled->Integral() << "    " << 1.0*PostTrigger_PseudoCaloMET_rescaled->Integral()/Nosel_PseudoCaloMET->Integral() << endl;
    cout << "After trigger + SF:  " << PostTrigger_PseudoCaloMET_SF_rescaled___Gluino->Integral() << "    " << 1.0*PostTrigger_PseudoCaloMET_SF_rescaled___Gluino->Integral()/Nosel_PseudoCaloMET->Integral() << endl;
    cout << "After PseudoMET>170: " << PostTrigger_PseudoCaloMET_SF_rescaled___Gluino_clone->Integral() << "    " << 1.0*PostTrigger_PseudoCaloMET_SF_rescaled___Gluino_clone->Integral()/Nosel_PseudoCaloMET->Integral() << endl;
    cout << "After HSCP sel:      " << METanalysis_Eta2p4_PseudoCaloMET->Integral() << "    " << 1.0*METanalysis_Eta2p4_PseudoCaloMET->Integral()/Nosel_PseudoCaloMET->Integral() << endl;


    return;
}


void CompareKinematics(const char *data = "JetMET2024F_V12p313.root",
                       const char *gluino1400 = "../output/Gluino_V19/Gluino_Run3_MET_madgraph_1400_V19p6.root",
                       const char *gluino2000 = "../output/Gluino_V19/Gluino_Run3_MET_madgraph_2000_V19p6.root",
                       const char *gluino2600 = "../output/Gluino_V19/Gluino_Run3_MET_madgraph_2600_V19p6.root") {

    TFile *ifile = new TFile(data, "READ");

    // HSCP-matched signal, one per mass slot; an empty path leaves the slot out
    struct Signal { const char *path; const char *key; const char *label; const char *tag; TH2F *h; };
    std::vector<Signal> signals;
    for (Signal s : std::vector<Signal>{{gluino1400, "gluino1400", "m#scale[0.7]{#tilde{g}}=1400 GeV", "1400", nullptr},
                                        {gluino2000, "gluino2000", "m#scale[0.7]{#tilde{g}}=2000 GeV", "2000", nullptr},
                                        {gluino2600, "gluino2600", "m#scale[0.7]{#tilde{g}}=2600 GeV", "2600", nullptr}}) {
        TFile *f = OpenOptional(s.path);
        if (!f) { std::cout << "CompareKinematics: " << s.key << " left out" << std::endl; continue; }
        s.h = (TH2F*)f->Get("trackPT_vs_trackPseudoTrackPT__HSCPmatched");
        if (!s.h) Abort(Form("CompareKinematics: no trackPT_vs_trackPseudoTrackPT__HSCPmatched in %s", s.path));
        signals.push_back(s);
    }

    TH2F *trackPT_vs_trackPseudoTrackPT = (TH2F*)ifile->Get("trackPT_vs_trackPseudoTrackPT");
    TH2F *trackETA_vs_trackPseudoTrackETA = (TH2F*)ifile->Get("trackETA_vs_trackPseudoTrackETA");
    TH2F *trackPHI_vs_trackPseudoTrackPHI = (TH2F*)ifile->Get("trackPHI_vs_trackPseudoTrackPHI");

    TH2F *trackPT_vs_trackPseudoTrackPT__PFmuon = (TH2F*)ifile->Get("trackPT_vs_trackPseudoTrackPT__PFmuon");
    TH2F *trackPT_vs_trackPseudoTrackPT__PFpion = (TH2F*)ifile->Get("trackPT_vs_trackPseudoTrackPT__PFpion");


    trackPT_vs_trackPseudoTrackPT->GetXaxis()->SetTitle("PF p_{T} [GeV/c]");
    trackPT_vs_trackPseudoTrackPT->GetYaxis()->SetTitle("Track p_{T} [GeV/c]");
    trackPT_vs_trackPseudoTrackPT->GetYaxis()->SetTitleOffset(1.2);
    trackPT_vs_trackPseudoTrackPT->GetZaxis()->SetTitle("Events");
    trackPT_vs_trackPseudoTrackPT->SetTitle("");
    trackPT_vs_trackPseudoTrackPT->GetYaxis()->SetTitleSize(0.06);
    trackPT_vs_trackPseudoTrackPT->GetXaxis()->SetTitleSize(0.06);
    trackPT_vs_trackPseudoTrackPT->GetXaxis()->SetTitleOffset(0.9);
    trackPT_vs_trackPseudoTrackPT->GetYaxis()->SetTitleOffset(1);
    trackPT_vs_trackPseudoTrackPT->GetXaxis()->SetLabelSize(0.05);
    trackPT_vs_trackPseudoTrackPT->GetYaxis()->SetLabelSize(0.05);
    trackPT_vs_trackPseudoTrackPT->GetZaxis()->SetTitleSize(0.06);
    trackPT_vs_trackPseudoTrackPT->GetZaxis()->SetTitleOffset(0.9);
    trackPT_vs_trackPseudoTrackPT->GetZaxis()->SetLabelSize(0.05);

    trackETA_vs_trackPseudoTrackETA->GetXaxis()->SetTitle("PF #eta");
    trackETA_vs_trackPseudoTrackETA->GetYaxis()->SetTitle("Track #eta");
    trackETA_vs_trackPseudoTrackETA->GetZaxis()->SetTitle("Events");
    trackETA_vs_trackPseudoTrackETA->SetTitle("");
    trackETA_vs_trackPseudoTrackETA->GetYaxis()->SetTitleSize(0.06);
    trackETA_vs_trackPseudoTrackETA->GetXaxis()->SetTitleSize(0.06);
    trackETA_vs_trackPseudoTrackETA->GetXaxis()->SetTitleOffset(0.9);
    trackETA_vs_trackPseudoTrackETA->GetYaxis()->SetTitleOffset(0.7);
    trackETA_vs_trackPseudoTrackETA->GetXaxis()->SetLabelSize(0.05);
    trackETA_vs_trackPseudoTrackETA->GetYaxis()->SetLabelSize(0.05);
    trackETA_vs_trackPseudoTrackETA->GetZaxis()->SetTitleSize(0.06);
    trackETA_vs_trackPseudoTrackETA->GetZaxis()->SetTitleOffset(0.9);
    trackETA_vs_trackPseudoTrackETA->GetZaxis()->SetLabelSize(0.05);

    trackPHI_vs_trackPseudoTrackPHI->GetXaxis()->SetTitle("PF #phi (rad)");
    trackPHI_vs_trackPseudoTrackPHI->GetYaxis()->SetTitle("Track #phi (rad)");
    trackPHI_vs_trackPseudoTrackPHI->GetZaxis()->SetTitle("Events");
    trackPHI_vs_trackPseudoTrackPHI->SetTitle("");
    trackPHI_vs_trackPseudoTrackPHI->GetYaxis()->SetTitleSize(0.06);
    trackPHI_vs_trackPseudoTrackPHI->GetXaxis()->SetTitleSize(0.06);
    trackPHI_vs_trackPseudoTrackPHI->GetXaxis()->SetTitleOffset(0.9);
    trackPHI_vs_trackPseudoTrackPHI->GetYaxis()->SetTitleOffset(0.7);
    trackPHI_vs_trackPseudoTrackPHI->GetXaxis()->SetLabelSize(0.05);
    trackPHI_vs_trackPseudoTrackPHI->GetYaxis()->SetLabelSize(0.05);
    trackPHI_vs_trackPseudoTrackPHI->GetZaxis()->SetTitleSize(0.06);
    trackPHI_vs_trackPseudoTrackPHI->GetZaxis()->SetTitleOffset(0.9);
    trackPHI_vs_trackPseudoTrackPHI->GetZaxis()->SetLabelSize(0.05);

    // --- Styling helper for the extra PT histograms (same style as trackPT_vs_trackPseudoTrackPT) ---
    auto stylePT = [](TH2F *h, double sizelabel=0.04) {
        h->GetXaxis()->SetTitle("PF p_{T} [GeV/c]");
        h->GetYaxis()->SetTitle("Track p_{T} [GeV/c]");
        h->GetZaxis()->SetTitle("Events");
        h->SetTitle("");
        h->GetYaxis()->SetTitleSize(0.06);
        h->GetXaxis()->SetTitleSize(0.06);
        h->GetXaxis()->SetTitleOffset(0.9);
        h->GetYaxis()->SetTitleOffset(1);
        h->GetXaxis()->SetLabelSize(sizelabel);
        h->GetYaxis()->SetLabelSize(sizelabel);
        h->GetZaxis()->SetTitleSize(0.06);
        h->GetZaxis()->SetTitleOffset(0.9);
        h->GetZaxis()->SetLabelSize(0.05);
    };

    stylePT(trackPT_vs_trackPseudoTrackPT__PFmuon);
    stylePT(trackPT_vs_trackPseudoTrackPT__PFpion);
    for (auto &s : signals) stylePT(s.h);

    TLatex *tex = new TLatex(0.62, 0.91, Label("lumi", "109 fb^{-1} (13.6 TeV)"));
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);

    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);

    TLatex *latexPF = new TLatex(0.62, 0.95, "#bf{PF muons}");
    latexPF->SetNDC();
    latexPF->SetTextFont(42);
    latexPF->SetTextSize(0.04);

    TLatex *latexglu = new TLatex(0.85, 0.91, "#bf{PF muons}");
    latexglu->SetNDC();
    latexglu->SetTextFont(42);
    latexglu->SetTextSize(0.04);
    

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    
    c1->SetLeftMargin(0.16); c1->SetBottomMargin(0.16);c1->SetRightMargin(0.16);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);
    trackPT_vs_trackPseudoTrackPT->Draw("COLZ");
    tex->Draw();
    c1->SetLogz();
    latex1->Draw();
    c1->SaveAs(OutPath("CompareKinematics_PT.pdf"));
    latex1->SetTitle("#it{Private work (CMS data)}");
    c1->Modified();
    c1->Update();
    c1->SaveAs(OutPath("CompareKinematics_PT_bis.pdf"));
    latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Work in progress}");


    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    c2->SetRightMargin(0.16);
    c2->SetLeftMargin(0.16); c2->SetBottomMargin(0.16);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);
    trackETA_vs_trackPseudoTrackETA->Draw("COLZ");
    tex->Draw();
    c2->SetLogz();
    latex1->Draw();
    c2->SaveAs(OutPath("CompareKinematics_ETA.pdf"));
    latex1->SetTitle("#it{Private work (CMS data)}");
    c2->Modified();
    c2->Update();
    c2->SaveAs(OutPath("CompareKinematics_ETA_bis.pdf"));
    latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Work in progress}");

    TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
    c3->SetRightMargin(0.16);
    c3->SetLeftMargin(0.16); c3->SetBottomMargin(0.16);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);
    trackPHI_vs_trackPseudoTrackPHI->Draw("COLZ");
    tex->Draw();
    c3->SetLogz();
    latex1->Draw();
    c3->SaveAs(OutPath("CompareKinematics_PHI.pdf"));
    latex1->SetTitle("#it{Private work (CMS data)}");
    c3->Modified();
    c3->Update();
    c3->SaveAs(OutPath("CompareKinematics_PHI_bis.pdf"));
    latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Work in progress}");


    TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);
    c4->SetRightMargin(0.16);
    c4->SetLeftMargin(0.16); c4->SetBottomMargin(0.16);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);
    trackPT_vs_trackPseudoTrackPT__PFmuon->Draw("COLZ");
    tex->Draw();
    c4->SetLogz();
    latex1->Draw();
    latexPF->Draw();
    c4->SaveAs(OutPath("CompareKinematics_PT_PFmuon.pdf"));
    latex1->SetTitle("#it{Private work (CMS data)}");
    c4->Modified();
    c4->Update();
    c4->SaveAs(OutPath("CompareKinematics_PT_PFmuon_bis.pdf"));
    latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Work in progress}");


    TCanvas *c5 = new TCanvas("c5", "c5", 800, 600);
    c5->SetRightMargin(0.16);
    c5->SetLeftMargin(0.16); c5->SetBottomMargin(0.16);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);
    trackPT_vs_trackPseudoTrackPT__PFpion->Draw("COLZ");
    tex->Draw();
    c5->SetLogz();
    latex1->Draw();
    latexPF->SetTitle("#bf{PF pions}");
    latexPF->Draw();
    c5->SaveAs(OutPath("CompareKinematics_PT_PFpion.pdf"));
    latex1->SetTitle("#it{Private work (CMS data)}");
    c5->Modified();
    c5->Update();
    c5->SaveAs(OutPath("CompareKinematics_PT_PFpion_bis.pdf"));
    latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Work in progress}");


    // --- HSCP matched: one canvas per mass slot given ---
    for (auto &s : signals) {
        TCanvas *c = new TCanvas(Form("c_HSCPmatched_%s", s.tag), Form("c_HSCPmatched_%s", s.tag), 800, 600);
        c->SetRightMargin(0.16);
        c->SetLeftMargin(0.16); c->SetBottomMargin(0.16);
        gStyle->SetOptStat(0);
        gStyle->SetPalette(kViridis);
        s.h->Draw("COLZ");
        c->SetLogz();
        latex1->Draw();
        latexglu->SetTitle(Form("#bf{HSCP matched (%s)}", Label(s.key, s.label)));
        latexglu->SetX(0.53);
        latexglu->Draw();
        c->SaveAs(OutPath("CompareKinematics_PT_HSCPmatched_%s.pdf", s.tag));
        latex1->SetTitle("#it{Private work (CMS simulation)}");
        c->Modified();
        c->Update();
        c->SaveAs(OutPath("CompareKinematics_PT_HSCPmatched_%s_bis.pdf", s.tag));
        latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    }


    return;
}

void LangausFitOnIh() {

    TFile *ofile = new TFile("../outputDisplay/LangausFitOnIh.root", "RECREATE");

    TFile *ifileWjets = new TFile("../output/Wjets2024_V14/WjetMuNu2024_V14p10_weighted.root", "READ");
    TFile *ifileJetMET = new TFile("../output/JetMET2024_V12/JetMET2024_V12p24.root", "READ");

    TH1F *Ih_wjet = (TH1F*)ifileWjets->Get("Nm1_event_Ih_StripOnly");
    TH1F *Ih_jetmet = (TH1F*)ifileJetMET->Get("Nm1_event_Ih_StripOnly");

    Ih_wjet->SetLineColor(kBlue-7);
    Ih_wjet->SetFillColorAlpha(kBlue-7, 0.5);
    Ih_jetmet->SetLineColor(kBlack);
    Ih_jetmet->SetMarkerColor(kBlack);
    Ih_jetmet->SetMarkerStyle(20);
    Ih_wjet->GetXaxis()->SetTitle("I_{h} (StripOnly)");
    Ih_wjet->GetYaxis()->SetTitle("Events");

    Ih_wjet->Scale(Ih_jetmet->Integral()/Ih_wjet->Integral());


    /////////////////////////////////////////////////////////

    // W+jets
    double peak_wjet  = Ih_wjet->GetBinCenter(Ih_wjet->GetMaximumBin());
    double xmin_wjet  = peak_wjet * 0.8;
    double xmax_wjet  = peak_wjet * 1.2;
    double area_wjet  = Ih_wjet->Integral() * Ih_wjet->GetBinWidth(1);

    TF1 *Langaus_wjet = new TF1("Langaus_wjet", langaufun, xmin_wjet, xmax_wjet, 4);
    Langaus_wjet->SetParNames("Width", "MPV", "Area", "GSigma");
    Langaus_wjet->SetParameters(0.05, peak_wjet, area_wjet, 0.05);
    Langaus_wjet->SetParLimits(0, 1e-4, 1.0);
    Langaus_wjet->SetParLimits(1, xmin_wjet, xmax_wjet);
    Langaus_wjet->SetParLimits(2, 0.0, 1e9);
    Langaus_wjet->SetParLimits(3, 1e-4, 1.0);
    Ih_wjet->Fit("Langaus_wjet", "RMQ");

    // JetMET
    double peak_jm   = Ih_jetmet->GetBinCenter(Ih_jetmet->GetMaximumBin());
    double xmin_jm   = peak_jm * 0.8;
    double xmax_jm   = peak_jm * 1.2;
    double area_jm   = Ih_jetmet->Integral() * Ih_jetmet->GetBinWidth(1);

    TF1 *Langaus_jetmet = new TF1("Langaus_jetmet", langaufun, xmin_jm, xmax_jm, 4);
    Langaus_jetmet->SetParNames("Width", "MPV", "Area", "GSigma");
    Langaus_jetmet->SetParameters(0.05, peak_jm, area_jm, 0.05);
    Langaus_jetmet->SetParLimits(0, 1e-4, 1.0);
    Langaus_jetmet->SetParLimits(1, xmin_jm, xmax_jm);
    Langaus_jetmet->SetParLimits(2, 0.0, 1e9);
    Langaus_jetmet->SetParLimits(3, 1e-4, 1.0);
    Ih_jetmet->Fit("Langaus_jetmet", "RMQ");

    Langaus_wjet->SetLineColor(kGreen);
    Langaus_jetmet->SetLineColor(kRed);


    // Drawing
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    gStyle->SetOptStat(0);
    Ih_wjet->Draw("HIST");
    Ih_jetmet->Draw("E1 same");
    Langaus_jetmet->Draw("same");
    Langaus_wjet->Draw("same");
    TLegend *legend = new TLegend(0.6, 0.7, 0.85, 0.9);
    legend->AddEntry(Ih_wjet, "W+jets", "f");
    legend->AddEntry(Ih_jetmet, "JetMET", "lep");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->Draw();
    c1->SetLogy();


    cout << "\t\t\t MPV" << endl;
    cout << "W+jets: " << Langaus_wjet->GetParameter("MPV") << " +/- " << Langaus_wjet->GetParError(Langaus_wjet->GetParNumber("MPV")) << endl;
    cout << "JetMET: " << Langaus_jetmet->GetParameter("MPV") << " +/- " << Langaus_jetmet->GetParError(Langaus_jetmet->GetParNumber("MPV")) << endl;

    // OUTPUT: data/MC = 3.2168/3.12692

    ofile->cd();
    c1->Write();
    ofile->Close();

    c1->SaveAs("../outputDisplay/LangausFitOnIh.pdf");


    return;
}

void ShowSignalEfficiency(const char *gluinoPattern = "../output/Gluino_V19/Gluino_Run3_MET_madgraph_%d_V19p3_weighted.root") {

    // one weighted file per mass: gluinoPattern with %d replaced by the mass; masses without a file are skipped
    std::map<int, TFile*> ifiles;
    for (int m : {1100, 1200, 1300, 1400, 1600, 1800, 2000, 2200, 2400, 2600}) {
        TString path = TString::Format(gluinoPattern, m);
        if (gSystem->AccessPathName(path)) {
            cout << "ShowSignalEfficiency: no file for mass " << m << " (" << path << "), skipped" << endl;
            continue;
        }
        ifiles[m] = new TFile(path, "READ");
    }

    std::vector <string> etaname = {"METanalysis_PseudoMETrescaled_Eta1",
                            "METanalysis_PseudoMETrescaled_Eta1_2p4",
                            "METanalysis_PseudoMETrescaled_Eta2p4"};

    TH1D* pEff_eta1 = new TH1D("eff_eta1", ";Gluino mass [GeV];Signal efficiency", 16, 1050, 2650);
    TH1D* pEff_eta1_2p4 = new TH1D("eff_eta1_2p4", ";Gluino mass [GeV];Signal efficiency", 16, 1050, 2650);
    TH1D* pEff_eta2p4 = new TH1D("eff_eta2p4", ";Gluino mass [GeV];Signal efficiency", 16, 1050, 2650);
    std::vector<TH1D*> pEffs = {pEff_eta1, pEff_eta1_2p4, pEff_eta2p4};
    for (int j = 0; j < etaname.size(); j++) {
        
        for (int i = 1100; i <= 2600; i += 100) {
            if (!ifiles.count(i)) continue;
            TFile *ifile = ifiles[i];

            TH1D *WithSelection = (TH1D*)ifile->Get(Form("%s_nHSCP", etaname[j].c_str()));
            TH1D *WithoutSelection = (TH1D*)ifile->Get("nHSCP");

            int bin = (i - 1050) / 100 + 1;
            pEffs[j]->SetBinContent(bin, WithSelection->Integral()/WithoutSelection->Integral());
            pEffs[j]->SetBinError(bin, sqrt(WithSelection->Integral())/WithoutSelection->Integral());
        }

    }

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    gStyle->SetOptStat(0);
    pEffs[0]->SetLineColor(kRed);
    pEffs[0]->SetMarkerColor(kRed);
    pEffs[0]->SetMarkerStyle(20);
    pEffs[0]->Draw("E1");

    pEffs[1]->SetLineColor(kGreen-3);
    pEffs[1]->SetMarkerColor(kGreen-3);
    pEffs[1]->SetMarkerStyle(21);
    pEffs[1]->Draw("E1 same");
    
    pEffs[2]->SetLineColor(kBlue-7);
    pEffs[2]->SetMarkerColor(kBlue-7);
    pEffs[2]->SetMarkerStyle(22);
    pEffs[2]->Draw("E1 same");
    

    pEffs[0]->SetTitle("");
    pEffs[0]->GetYaxis()->SetRangeUser(0.0, 1.0);
    TLegend *legend = new TLegend(0.6, 0.7, 0.85, 0.9);
    legend->AddEntry(pEffs[0], "|#eta| < 1.0", "lep");
    legend->AddEntry(pEffs[1], "1.0 < |#eta| < 2.4", "lep");
    legend->AddEntry(pEffs[2], "|#eta| < 2.4", "lep");
    legend->Draw();

    c1->SaveAs(OutPath("SignalEfficiency_vs_GluinoMass.pdf"));


    return;
}


void FpixSlices(const char *wjets = "../output/Wjets2024_V14/WjetMuNu2024_V14p11_weighted.root",
                const char *ttbar = "../output/TTbar2024_V15/TTbar2024_V15p5_weighted.root",
                const char *qcd = "../output/QCD2024_V16/QCD2024_mu_V16p1_weighted.root",
                const char *gluino = "../output/Gluino_V19/Gluino_Run3_MET_madgraph_2000_V19p3.root") {

    // backgrounds: an empty path leaves the sample out of the MC sum
    TFile *fileWjets = OpenOptional(wjets);
    TFile *fileTTbar = OpenOptional(ttbar);
    TFile *fileQCD = OpenOptional(qcd);
    TFile *fileSignal = new TFile(gluino, "READ");
    
    TH1F *FpixGluino_Eta1 = (TH1F*)fileSignal->Get("METanalysis_PseudoMETrescaled_Eta1_Fpix");
    TH1F *FpixGluino_Eta1_2p4 = (TH1F*)fileSignal->Get("METanalysis_PseudoMETrescaled_Eta1_2p4_Fpix");
    TH1F *FpixGluino_Eta2p4 = (TH1F*)fileSignal->Get("METanalysis_PseudoMETrescaled_Eta2p4_Fpix");

    // per background given: Fpix for |eta| < 1, 1 < |eta| < 2.4, |eta| < 2.4
    std::vector<std::vector<TH1F*>> mc;
    if (fileQCD) mc.push_back({(TH1F*)fileQCD->Get("METanalysis_Eta1_Fpix"),
                               (TH1F*)fileQCD->Get("METanalysis_Eta1_2p4_Fpix"),
                               (TH1F*)fileQCD->Get("METanalysis_Eta2p4_Fpix")});
    if (fileTTbar) mc.push_back({(TH1F*)fileTTbar->Get("METanalysis_Eta1_Fpix"),
                                 (TH1F*)fileTTbar->Get("METanalysis_Eta1_2p4_Fpix"),
                                 (TH1F*)fileTTbar->Get("METanalysis_Eta2p4_Fpix")});
    if (fileWjets) mc.push_back({(TH1F*)fileWjets->Get("METanalysis_PseudoMETrescaled_Eta1_Fpix"),
                                 (TH1F*)fileWjets->Get("METanalysis_PseudoMETrescaled_Eta1_2p4_Fpix"),
                                 (TH1F*)fileWjets->Get("METanalysis_PseudoMETrescaled_Eta2p4_Fpix")});
    if (mc.empty()) Abort("FpixSlices: give at least one background");
    for (auto &hs : mc)
        for (TH1F *h : hs)
            if (!h) Abort("FpixSlices: a background Fpix histogram is missing");

    // MC = sum of the backgrounds given (QCD + TTbar + Wjets)
    TH1F *FpixQCD_Eta1 = (TH1F*)mc[0][0]->Clone("FpixMC_Eta1");
    TH1F *FpixQCD_Eta1_2p4 = (TH1F*)mc[0][1]->Clone("FpixMC_Eta1_2p4");
    TH1F *FpixQCD_Eta2p4 = (TH1F*)mc[0][2]->Clone("FpixMC_Eta2p4");
    for (size_t k = 1; k < mc.size(); k++) {
        FpixQCD_Eta1->Add(mc[k][0]);
        FpixQCD_Eta1_2p4->Add(mc[k][1]);
        FpixQCD_Eta2p4->Add(mc[k][2]);
    }
    

    for (int i = 1; i <= 10; i++) {
        double prob = i * 0.1;
        double qGluino_Eta2p4, qGluino_Eta1, qGluino_Eta1_2p4, qMC_eta2p4, qMC_eta1, qMC_eta1_2p4;

        FpixGluino_Eta2p4->GetQuantiles(1, &qGluino_Eta2p4, &prob);
        FpixQCD_Eta2p4   ->GetQuantiles(1, &qMC_eta2p4,    &prob);
        FpixGluino_Eta1->GetQuantiles(1, &qGluino_Eta1, &prob);
        FpixQCD_Eta1   ->GetQuantiles(1, &qMC_eta1,    &prob);
        FpixGluino_Eta1_2p4->GetQuantiles(1, &qGluino_Eta1_2p4, &prob);
        FpixQCD_Eta1_2p4   ->GetQuantiles(1, &qMC_eta1_2p4,    &prob);

        cout << "      Quantile " << prob << ":\n";
        cout << "Gluino |#eta| < 2.4: " << qGluino_Eta2p4 << "\n";
        cout << "Gluino |#eta| < 1.0: " << qGluino_Eta1 << "\n";
        cout << "Gluino 1.0 < |#eta| < 2.4: " << qGluino_Eta1_2p4 << "\n";
        cout << "MC (QCD+TTbar+Wjets) |#eta| < 2.4: " << qMC_eta2p4 << "\n";
        cout << "MC (QCD+TTbar+Wjets) |#eta| < 1.0: " << qMC_eta1 << "\n";
        cout << "MC (QCD+TTbar+Wjets) 1.0 < |#eta| < 2.4: " << qMC_eta1_2p4 << "\n";
    }

    
    
    return;
}


void DefineIhCut() {

    TFile *ifileData = new TFile("../output/JetMET2024_V12/JetMET2024_V12p24.root", "READ");

    TH1F *Ih = (TH1F*)ifileData->Get("Nm1_Ih_StripOnly");

    double C_2024 = 2.9784;


    double totalIntegral = Ih->Integral(C_2024, Ih->GetNbinsX() + 1);
    double targetIntegral = totalIntegral * 0.01;
    double cumulativeIntegral = 0.0;
    double x_b = 0.0;
    for (int i = Ih->GetNbinsX(); i >= 1; i--) {
        cumulativeIntegral += Ih->GetBinContent(i);
        if (cumulativeIntegral >= targetIntegral) {
            x_b = Ih->GetBinLowEdge(i);
            break;
        }
    }

    Ih->SetLineColor(kBlue);
    
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    gStyle->SetOptStat(0);
    Ih->Draw("HIST");
    TLine *line = new TLine(x_b, 0, x_b, Ih->GetMaximum());
    line->SetLineColor(kRed);
    TLine *line2 = new TLine(C_2024, 0, C_2024, Ih->GetMaximum());
    line2->SetLineColor(kBlack);
    line->Draw("same");
    line2->Draw("same");
    c1->SetLogy();

    c1->SaveAs("../outputDisplay/DefineIhCut.pdf");

    // check


    cout << endl;
    cout << "   xb: " << x_b << endl; 
    cout << "Integral above C: " << totalIntegral << endl;
    cout << "cumulativeIntegral: " << cumulativeIntegral << endl;
    cout << "ratio: " << cumulativeIntegral / totalIntegral << endl;

    return;
}


void SignalAcceptance_EtaSlice(bool isHSCPcharged = false) {

    const std::string version = isHSCPcharged ? "19p4" : "19p8";

    const int nSlices = 24;  // |eta| de 0.1 à 2.4 par pas de 0.1

    const char *nameNosel = isHSCPcharged ? "Nosel_GenHSCPcharged_Eta" : "Nosel_eta";
    const char *nameSel   = isHSCPcharged ? "HSCPsel_GenHSCPcharged_Eta"
                                          : "METanalysis_TestPUppiMETCut_Eta2p4_eta";

    TFile *ofile = new TFile("../outputDisplay/SignalAcceptance_EtaSlice.root", "RECREATE");

    // --- lambda : construit l'histo d'acceptance a partir d'un fichier ---
    auto makeAcc = [&](TFile *f, const char *hname) -> TH1D* {
        if (!f || f->IsZombie()) { std::cerr << "Bad file for " << hname << std::endl; return nullptr; }
        TH1D *h_nosel = (TH1D*)f->Get(nameNosel);
        TH1D *h_sel   = (TH1D*)f->Get(nameSel);
        if (!h_nosel || !h_sel) { std::cerr << "Missing histos in " << f->GetName() << std::endl; return nullptr; }

        ofile->cd();
        TH1D *h_acc = new TH1D(hname, ";|#eta| < x;Acceptance", nSlices, 0.05, 2.45);
        h_acc->SetDirectory(ofile);

        for (int s = 1; s <= nSlices; s++) {
            double etaMax = 0.1 * s;
            int binLow  = h_sel->FindBin(-etaMax + 1e-6);
            int binHigh = h_sel->FindBin( etaMax - 1e-6);

            double errSel = 0., errNosel = 0.;
            double nSel   = h_sel->IntegralAndError(binLow, binHigh, errSel);
            double nNosel = h_nosel->IntegralAndError(0, h_nosel->GetNbinsX()+1, errNosel);

            double acc = (nNosel > 0) ? nSel / nNosel : 0.;
            double err = 0.;
            if (nSel > 0 && nNosel > 0)
                err = acc * std::sqrt(std::pow(errSel/nSel, 2) + std::pow(errNosel/nNosel, 2));

            h_acc->SetBinContent(s, acc);
            h_acc->SetBinError(s, err);
        }
        return h_acc;
    };

    // ================= Signal (3 masses) =================
    std::vector<int> masses = {2000, 2400, 2600};
    std::vector<TH1D*> h_accs;

    for (int mass : masses) {
        std::string fname = "../output/Gluino_V19/Gluino_Run3_MET_madgraph_"
                          + std::to_string(mass) + "_V" + version + "_weighted.root";
        TFile *f = new TFile(fname.c_str(), "READ");
        TH1D *h_acc = makeAcc(f, Form("Acceptance_EtaSlice_M%d", mass));
        if (h_acc) {
            std::cout << mass << " nosel integral = "
                      << ((TH1D*)f->Get(nameNosel))->Integral(0, ((TH1D*)f->Get(nameNosel))->GetNbinsX()) << std::endl;
            h_acc->Write();
            h_accs.push_back(h_acc);
        }
    }

    // ================= Data =================
    TFile *fData = new TFile("../output/JetMET2024_V12/JetMET2024_V12p32.root", "READ");
    TH1D *h_accData = makeAcc(fData, "Acceptance_EtaSlice_Data");
    if (h_accData) h_accData->Write();

    // --- Canvas récapitulatif ---
    TCanvas *c = new TCanvas("c_SignalAcceptance_EtaSlice", "Signal acceptance vs |#eta| cut", 800, 600);
    c->SetGrid();
    c->SetLeftMargin(0.16); c->SetBottomMargin(0.16);

    TLegend *leg = new TLegend(0.2, 0.65, 0.45, 0.89);
    leg->SetBorderSize(0);

    int colors[3] = {kOrange+8, kViolet+1, kGreen-3};
    int marker[3] = {22, 21, 23};

    double ymax = 0.;
    for (auto h : h_accs) ymax = std::max(ymax, h->GetMaximum());
    if (h_accData) ymax = std::max(ymax, h_accData->GetMaximum());

    auto styleAxes = [](TH1D *h) {
        h->SetStats(0);
        h->SetTitle("");
        h->SetLineWidth(2);
        h->GetYaxis()->SetTitleSize(0.06);
        h->GetXaxis()->SetTitleSize(0.06);
        h->GetXaxis()->SetTitleOffset(0.9);
        h->GetYaxis()->SetTitleOffset(1);
        h->GetXaxis()->SetLabelSize(0.05);
        h->GetYaxis()->SetLabelSize(0.05);
    };

    bool first = true;
    for (size_t i = 0; i < h_accs.size(); i++) {
        TH1D *h = h_accs[i];
        styleAxes(h);
        h->SetLineColor(colors[i]);
        h->SetMarkerColor(colors[i]);
        h->SetMarkerStyle(marker[i]);
        h->Scale(0.87);
        h->GetYaxis()->SetRangeUser(0., ymax);
        h->Draw(first ? "PE" : "PE SAME");
        first = false;
        leg->AddEntry(h, Form("M = %d GeV", masses[i]), "lp");
    }

    // data en noir, par-dessus
    if (h_accData) {
        styleAxes(h_accData);
        h_accData->SetLineColor(kBlack);
        h_accData->SetMarkerColor(kBlack);
        h_accData->SetMarkerStyle(20);
        h_accData->GetYaxis()->SetRangeUser(0., ymax * 1.2);
        h_accData->Scale(124);
        h_accData->Draw(first ? "PE" : "PE SAME");
        leg->AddEntry(h_accData, "Data (x400)", "lp");
    }

    leg->Draw();

    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);
    latex1->Draw();

    TLatex *latex2 = new TLatex(0.68, 0.91, "109 fb^{-1} (13.6 TeV)");
    latex2->SetNDC();
    latex2->SetTextFont(42);
    latex2->SetTextSize(0.04);
    latex2->Draw();

    c->SaveAs((isHSCPcharged)? "../outputDisplay/SignalAcceptance_EtaSlice_HSCPcharged.pdf" : "../outputDisplay/SignalAcceptance_EtaSlice.pdf");

    ofile->cd();
    c->Write();

    latex1->SetTitle("#it{Private work (CMS simulation/data)}");
    c->Modified();
    c->Update();
    c->SaveAs((isHSCPcharged)? "../outputDisplay/SignalAcceptance_EtaSlice_HSCPcharged_bis.pdf" : "../outputDisplay/SignalAcceptance_EtaSlice_bis.pdf");

    ofile->Close();

    return;
}


void Corr_Ih_1oP(std::string eta) {

    TFile *ifile_data = new TFile("../output/JetMET2024_V12/JetMET2024_V12p27.root", "READ");

    TH2F *h_10000oP_vs_Ih__nosel = (TH2F*)ifile_data->Get("Nosel_10000oP_vs_Ih");
    TH2F *h_10000oP_vs_Ih__sel   = (TH2F*)ifile_data->Get(Form("METanalysis_PseudoMETrescaled_%s_10000oP_vs_Ih", eta.c_str()));

    // make the profil of each and then draw
    TProfile *p_nsel = h_10000oP_vs_Ih__nosel->ProfileX();
    TProfile *p_sel = h_10000oP_vs_Ih__sel->ProfileX();

    h_10000oP_vs_Ih__nosel->GetXaxis()->SetTitle("10^{4}/p [GeV^{-1}]");
    h_10000oP_vs_Ih__nosel->GetYaxis()->SetTitle("I_{h} [MeV/cm]");
    h_10000oP_vs_Ih__nosel->GetZaxis()->SetTitle("Events");
    h_10000oP_vs_Ih__nosel->SetTitle("");
    h_10000oP_vs_Ih__nosel->GetYaxis()->SetTitleSize(0.06);
    h_10000oP_vs_Ih__nosel->GetXaxis()->SetTitleSize(0.06);
    h_10000oP_vs_Ih__nosel->GetXaxis()->SetTitleOffset(0.9);
    h_10000oP_vs_Ih__nosel->GetYaxis()->SetTitleOffset(0.9);
    h_10000oP_vs_Ih__nosel->GetXaxis()->SetLabelSize(0.05);
    h_10000oP_vs_Ih__nosel->GetYaxis()->SetLabelSize(0.05);
    h_10000oP_vs_Ih__nosel->GetZaxis()->SetTitleSize(0.06);
    h_10000oP_vs_Ih__nosel->GetZaxis()->SetTitleOffset(0.9);
    h_10000oP_vs_Ih__nosel->GetZaxis()->SetLabelSize(0.05);

    h_10000oP_vs_Ih__sel->GetXaxis()->SetTitle("10^{4}/p [GeV^{-1}]");
    h_10000oP_vs_Ih__sel->GetYaxis()->SetTitle("I_{h} [MeV/cm]");
    h_10000oP_vs_Ih__sel->GetZaxis()->SetTitle("Events");
    h_10000oP_vs_Ih__sel->SetTitle("");
    h_10000oP_vs_Ih__sel->GetYaxis()->SetTitleSize(0.06);
    h_10000oP_vs_Ih__sel->GetXaxis()->SetTitleSize(0.06);
    h_10000oP_vs_Ih__sel->GetXaxis()->SetTitleOffset(0.9);
    h_10000oP_vs_Ih__sel->GetYaxis()->SetTitleOffset(0.9);
    h_10000oP_vs_Ih__sel->GetXaxis()->SetLabelSize(0.05);
    h_10000oP_vs_Ih__sel->GetYaxis()->SetLabelSize(0.05);
    h_10000oP_vs_Ih__sel->GetZaxis()->SetTitleSize(0.06);
    h_10000oP_vs_Ih__sel->GetZaxis()->SetTitleOffset(0.9);
    h_10000oP_vs_Ih__sel->GetZaxis()->SetLabelSize(0.05);

    p_nsel->SetLineColor(kRed);
    p_nsel->SetLineWidth(2);
    p_sel->SetLineColor(kRed);
    p_sel->SetLineWidth(2);

    // --- fit linéaire (non dessiné) sur les profils ---
    TF1 *fit_nsel = new TF1("fit_nsel", "[0]+[1]*x", 0, 210);
    TF1 *fit_sel  = new TF1("fit_sel",  "[0]+[1]*x", 0, 210);

    p_nsel->Fit(fit_nsel, "RQ0");
    p_sel->Fit(fit_sel,  "RQ0");

    cout << "Nosel : I_h = " << fit_nsel->GetParameter(0) << " + " << fit_nsel->GetParameter(1) << " * (10^4/p)" << endl;
    cout << "Sel   : I_h = " << fit_sel->GetParameter(0)  << " + " << fit_sel->GetParameter(1)  << " * (10^4/p)" << endl;

    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);

    TLatex *tex = new TLatex(0.62, 0.91, "105.8 fb^{-1} (13.6 TeV)");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);

    TCanvas *c_nosel = new TCanvas("c_nosel", "c_nosel", 800, 600);
    c_nosel->SetGrid();
    c_nosel->SetLeftMargin(0.16); c_nosel->SetBottomMargin(0.16); c_nosel->SetRightMargin(0.19);
    h_10000oP_vs_Ih__nosel->Draw("COLZ");
    h_10000oP_vs_Ih__nosel->GetYaxis()->SetRangeUser(2.5, 6);
    h_10000oP_vs_Ih__nosel->GetXaxis()->SetRangeUser(0, 210);
    gStyle->SetOptStat(0);
    p_nsel->Draw("same");
    latex1->Draw();
    tex->Draw();
    c_nosel->SaveAs(Form("../outputDisplay/Corr_Ih_1oP__%s_nosel.pdf", eta.c_str()));
    latex1->SetTitle("#it{Private work (CMS data)}");
    c_nosel->Modified();
    c_nosel->Update();
    c_nosel->SaveAs(Form("../outputDisplay/Corr_Ih_1oP__%s_nosel_bis.pdf", eta.c_str()));

    TCanvas *c_sel = new TCanvas("c_sel", "c_sel", 800, 600);
    c_sel->SetGrid();
    c_sel->SetLeftMargin(0.16); c_sel->SetBottomMargin(0.16); c_sel->SetRightMargin(0.16);
    h_10000oP_vs_Ih__sel->Draw("COLZ");
    h_10000oP_vs_Ih__sel->GetYaxis()->SetRangeUser(2.5, 6);
    h_10000oP_vs_Ih__sel->GetXaxis()->SetRangeUser(0, 210);
    gStyle->SetOptStat(0);
    p_sel->Draw("same");
    latex1->Draw();
    tex->Draw();
    c_sel->SaveAs(Form("../outputDisplay/Corr_Ih_1oP__%s_sel.pdf", eta.c_str()));
    latex1->SetTitle("#it{Private work (CMS data)}");
    c_sel->Modified();
    c_sel->Update();
    c_sel->SaveAs(Form("../outputDisplay/Corr_Ih_1oP__%s_sel_bis.pdf", eta.c_str()));

    return;
}


void Acceptance_EtaSlice_DataMC(const char *data = "../output/JetMET2024_V12/JetMET2024_V12p27.root",
                                const char *wjets = "../output/Wjets2024_V14/WjetMuNu2024_V14p8_weighted.root",
                                const char *ttbar = "../output/TTbar2024_V15/TTbar2024_V15p6_weighted.root",
                                const char *qcd = "../output/QCD2024_V16/QCD2024_mu_V16p1_weighted.root") {

    // data first, then the backgrounds given (an empty path leaves the sample out)
    const char *paths[] = {data, wjets, ttbar, qcd};
    const char *legKeys[] = {"data", "wjets", "ttbar", "qcd"};
    const char *allLabels[] = {"Data", "W+jets", "t#bar{t}", "QCD"};
    std::vector<TFile*> files;
    std::vector<TString> labels;
    std::vector<const char*> keys;
    for (int k = 0; k < 4; k++) {
        TFile *f = (k == 0) ? new TFile(paths[k], "READ") : OpenOptional(paths[k]);
        if (!f) { std::cout << "Acceptance_EtaSlice_DataMC: " << legKeys[k] << " left out" << std::endl; continue; }
        files.push_back(f);
        labels.push_back(allLabels[k]);
        keys.push_back(legKeys[k]);
    }

    std::vector<TH1D*> h_accs;
    std::vector<int> accIndex;       // index in files/labels/keys of each entry of h_accs

    TFile *ofile = new TFile(OutPath("Acceptance_EtaSlice_DataMC.root"), "RECREATE");

    const int nSlices = 24;  // |eta| de 0.1 à 2.4 par pas de 0.1

    for (int i = 0; i < (int)files.size(); i++) {
        TH1D *h_eta_nosel = (TH1D*)files[i]->Get("Nosel_eta");
        TH1D *h_eta_sel   = (TH1D*)files[i]->Get((labels[i]=="Data" ? "METanalysis_PseudoMETnotRescaled_Eta2p4_eta" : "METanalysis_Eta2p4_eta"));

        if (!h_eta_nosel || !h_eta_sel) {
            std::cout << "Histogrammes manquants dans " << files[i]->GetName() << std::endl;
            continue;
        }

        ofile->cd();
        TH1D *h_acc = new TH1D(Form("Acceptance_EtaSlice_%s", labels[i].Data()),
                               Form("Acceptance vs |#eta| cut (%s);|#eta| < x;Acceptance", labels[i].Data()),
                               nSlices, 0.05, 2.45);  // bins centrés sur 0.1, 0.2, ..., 2.4

        // Dénominateur : intégrale totale (underflow + overflow inclus), calculé une fois
        double errNosel = 0.;
        double nNosel = h_eta_nosel->IntegralAndError(0, h_eta_nosel->GetNbinsX() + 1, errNosel);

        for (int s = 1; s <= nSlices; s++) {
            double etaMax = 0.1 * s;  // 0.1, 0.2, ..., 2.4

            int binLow_sel  = h_eta_sel->FindBin(-etaMax + 1e-6);
            int binHigh_sel = h_eta_sel->FindBin( etaMax - 1e-6);

            double errSel = 0.;
            double nSel = h_eta_sel->IntegralAndError(binLow_sel, binHigh_sel, errSel);

            double acc = (nNosel > 0) ? nSel / nNosel : 0.;
            double err = 0.;
            if (nSel > 0 && nNosel > 0)
                err = acc * std::sqrt(std::pow(errSel / nSel, 2) + std::pow(errNosel / nNosel, 2));

            h_acc->SetBinContent(s, acc);
            h_acc->SetBinError(s, err);
        }

        std::cout << labels[i] << " : Nosel integral = " << nNosel << std::endl;

        h_acc->Write();
        h_accs.push_back(h_acc);
        accIndex.push_back(i);
    }

    // --- Canvas récapitulatif ---
    TCanvas *c = new TCanvas("c_Acceptance_EtaSlice_DataMC", "Acceptance vs |#eta| cut", 800, 600);
    c->SetGrid();
    c->SetLeftMargin(0.16); c->SetBottomMargin(0.16);

    TLegend *leg = new TLegend(0.2, 0.7, 0.4, 0.89);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.05);

    // Data en noir, MC en couleurs
    int colors[4]  = {kBlack, kRed+1, kAzure+1, kGreen+2};
    int markers[4] = {20, 21, 22, 23};

    double ymax = 0.;
    for (auto h : h_accs) ymax = std::max(ymax, h->GetMaximum());

    for (size_t i = 0; i < h_accs.size(); i++) {
        TH1D *h = h_accs[i];

        h->SetLineColor(colors[i % 4]);
        h->SetMarkerColor(colors[i % 4]);
        h->SetMarkerStyle(markers[i % 4]);
        h->SetMarkerSize(0.8);
        h->SetLineWidth(2);
        h->SetStats(0);
        h->GetYaxis()->SetRangeUser(0., ymax * 1.2);
        h->SetTitle("");
        h->GetYaxis()->SetTitleSize(0.06);
        h->GetXaxis()->SetTitleSize(0.06);
        h->GetXaxis()->SetTitleOffset(0.9);
        h->GetYaxis()->SetTitleOffset(1.2);
        h->GetXaxis()->SetLabelSize(0.05);
        h->GetYaxis()->SetLabelSize(0.05);

        h->Draw(i == 0 ? "PE" : "PE SAME");
        leg->AddEntry(h, Label(keys[accIndex[i]], labels[accIndex[i]]), "lp");
    }
    leg->Draw();

    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{ Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);
    latex1->Draw();

    c->SaveAs(OutPath("Acceptance_EtaSlice_DataMC.pdf"));

    ofile->cd();
    c->Write();

    // Version "Private work"
    latex1->SetText(0.16, 0.91, "#it{Private work (CMS data)}");
    c->Modified();
    c->Update();
    c->SaveAs(OutPath("Acceptance_EtaSlice_DataMC_bis.pdf"));
    c->SaveAs(OutPath("Acceptance_EtaSlice_DataMC_bis.C"));

    ofile->Close();

    return;
}


//------------------------------------------------------------------
// Trigger efficiencies
//------------------------------------------------------------------s

void TriggerEffCalib__Signal(const char *labelSIGNAL, const char *inputfileSIGNAL, const char *ofilename) {

    gErrorIgnoreLevel = kError;

    TFile *ifileSIGNAL = new TFile(inputfileSIGNAL, "READ");

    if (!ifileSIGNAL || ifileSIGNAL->IsZombie()) {
        std::cerr << "Error: Could not open input file " << inputfileSIGNAL << std::endl;
        return;
    }

    TFile *ofile = new TFile(Form("../outputDisplay/TriggEff/TriggerEffCalib__Signal_%s.root", ofilename), "RECREATE");

    // ------------------------------------------------------------------
    // configuration
    // ------------------------------------------------------------------
    struct Observable {
        std::string label;                 // which selection label to read (nominal or rescaled)
        std::string branch;                // histo suffix in the file
        std::string tag;                   // short tag used in canvas/PDF names
        std::string xtitle;                // x axis title
        double cut;                        // lower threshold [GeV] for the "after-cut" efficiency (<0 = none)
        std::vector<std::string> triggers; // numerator trigger tags
    };

    std::vector<std::string> baseTriggers = {
        "HLT_PFMET120_PFMHT120_IDTight",
        "HLT_PFHT500_PFMET100_PFMHT100_IDTight",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",
        "HLT_MET105_IsoTrk50",
        "orMETtrg"
    };

    std::vector<std::string> orMET3a4 = {
        "orMET3a4trg1", "orMET3a4trg2", "orMET3a4trg3", "orMET3a4trg4"
    };

    // base triggers + the 4 orMET3a4trg
    std::vector<std::string> fullTriggers = baseTriggers;
    fullTriggers.insert(fullTriggers.end(), orMET3a4.begin(), orMET3a4.end());

    // marker sizes used in the original macro, indexed by base trigger
    std::map<std::string, double> markerSize = {
        { "HLT_PFMET120_PFMHT120_IDTight",                 0.045 },
        { "HLT_PFHT500_PFMET100_PFMHT100_IDTight",         0.04  },
        { "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",  0.032 },
        { "HLT_MET105_IsoTrk50",                           0.05  },
        { "orMETtrg",                                      0.06  }
    };

    std::string lbl     = labelSIGNAL;

    std::vector<Observable> observables = {
        // label,  branch,          tag,             xtitle,                   cut,   triggers
        { lbl,   "PseudoCaloMET", "PseudoCaloMET", "Pseudo MET [GeV]",       250.,  fullTriggers },
        { lbl,   "PUppiMET",      "PUppiMET",      "PUppi MET [GeV]",        150.,  fullTriggers },
        { lbl,   "PUppiMETNoMu",  "PUppiMETNoMu",  "PUppi MET (NoMu) [GeV]", 150.,  fullTriggers },
        { lbl,   "PFtrackPT",     "PFtrackPT",     "PF track p_{T} [GeV]",   -1.,   { "HLT_MET105_IsoTrk50" } }
    };

    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);

    // axis upper edge: PFtrackPT goes to 3000, MET observables to 2500
    // (drawing range kept at original-style values; tweak if needed)

    // ------------------------------------------------------------------
    // integrated-efficiency table (LaTeX)
    // ------------------------------------------------------------------
    std::ostringstream texBody;  // accumulates the table rows

    // binomial efficiency + error from integrated counts.
    // counts are summed over bins with low edge >= xmin (xmin < 0 -> full range,
    // under/overflow included).
    auto integratedEff = [](TH1F *num, TH1F *den, double xmin,
                            double &eff, double &err) {
        int blo, bhi;
        if (xmin < 0) {
            blo = 0;                       // include underflow
            bhi = den->GetNbinsX() + 1;    // include overflow
        } else {
            blo = den->GetXaxis()->FindBin(xmin);
            bhi = den->GetNbinsX() + 1;    // up to overflow
        }
        double N = den->Integral(blo, bhi);
        double k = num->Integral(blo, bhi);
        if (N > 0) {
            eff = k / N;
            double var = eff * (1.0 - eff) / N;   // gaussian binomial approx
            err = (var > 0) ? std::sqrt(var) : 0.0;
        } else {
            eff = 0.0;
            err = 0.0;
        }
    };

    // ------------------------------------------------------------------
    // loop over observables / triggers
    // ------------------------------------------------------------------
    for (const auto &obs : observables) {

        // inclusive denominator
        TH1F *den = (TH1F*)ifileSIGNAL->Get(Form("%s_%s", obs.label.c_str(), obs.branch.c_str()));

        if (!den) {
            std::cerr << "Warning: missing denominator " << obs.label << "_" << obs.branch
                      << ", skipping observable." << std::endl;
            continue;
        }

        // drawing x-range upper edge depends on the observable
        double xUp = (obs.branch == "PFtrackPT") ? 3000 : 1200;

        for (const std::string &trg : obs.triggers) {

            std::string numName = Form("%s_if___%s___%s",
                                       obs.label.c_str(), trg.c_str(), obs.branch.c_str());

            TH1F *num = (TH1F*)ifileSIGNAL->Get(numName.c_str());

            if (!num) {
                std::cerr << "Warning: missing numerator " << numName
                          << ", skipping." << std::endl;
                continue;
            }

            if (obs.branch == "PFtrackPT") { num->Rebin(3); den->Rebin(3); }

            // ---- efficiency = numerator / denominator ----
            std::string effName = Form("eff_%s_%s_SIGNAL", trg.c_str(), obs.tag.c_str());
            TH1F *eff = (TH1F*)num->Clone(effName.c_str());
            eff->SetDirectory(0);
            eff->Divide(den);

            if (eff->GetEntries() == 0 || eff->Integral() == 0) {
                std::cerr << "Warning: empty efficiency histo " << effName << ", skipping draw." << std::endl;
                continue;
            }

            // ---- integrated efficiency for the LaTeX table ----
            // full range
            double effInt = 0., effErr = 0.;
            integratedEff(num, den, -1., effInt, effErr);

            // after a lower cut on the observable (if defined for this obs)
            double effCut = -1., effCutErr = 0.;
            if (obs.cut >= 0)
                integratedEff(num, den, obs.cut, effCut, effCutErr);

            // escape underscores for LaTeX
            std::string obsTex = obs.tag; std::string trgTex = trg;
            for (auto *s : { &obsTex, &trgTex })
                for (size_t p = 0; (p = s->find('_', p)) != std::string::npos; p += 2)
                    s->replace(p, 1, "\\_");

            // cut column: value if a cut is defined, dash otherwise
            std::string cutCell;
            if (obs.cut >= 0)
                cutCell = Form("%.4f $\\pm$ %.4f", effCut, effCutErr);
            else
                cutCell = "--";

            texBody << "\\texttt{" << trgTex << "} & \\texttt{" << obsTex << "} & " << cutCell << " \\\\\n";

            // ---- drawing ----
            double msize = markerSize.count(trg) ? markerSize[trg] : 0.045;
            std::string cName  = Form("c_%s___%s", trg.c_str(), obs.tag.c_str());
            std::string ytitle = (trg == "orMETtrg") ? "eff orMETtrg"
                                                      : Form("eff. %s", trg.c_str());

            TCanvas *c = DrawCanvas(eff, cName.c_str(), obs.xtitle.c_str(), ytitle.c_str(),
                                    "E1", msize, 0, xUp, kBlack, 0, 1, false);
            c->cd();
            latex1->Draw();

            // ---- saving ----
            c->SaveAs(Form("../outputDisplay/TriggEff/c_%s___%s__SIGNAL__%s.pdf",
                           trg.c_str(), obs.tag.c_str(), ofilename));
            c->SaveAs(Form("../outputDisplay/TriggEff/c_%s___%s__SIGNAL__%s.C",
                           trg.c_str(), obs.tag.c_str(), ofilename));

            ofile->cd();
            c->Write();
        }
    }

    // ------------------------------------------------------------------
    // write LaTeX table
    // ------------------------------------------------------------------
    std::ofstream tex(Form("../outputDisplay/TriggEff/TriggerEff_table_%s.txt", ofilename));
    tex << "\\begin{table}[htbp]\n  \\centering\n"
        << "  \\caption{Integrated trigger efficiency (signal). "
        << "The last column gives the efficiency after a lower cut on the observable "
        << "(Pseudo MET $>$ 250 GeV, PUppi MET $>$ 150 GeV).}\n"
        << "  \\begin{tabular}{llc}\n    \\hline\n"
        << "    Trigger & Observable & $\\varepsilon$ (after cut) \\\\\n    \\hline\n"
        << texBody.str()
        << "    \\hline\n  \\end{tabular}\n\\end{table}\n";
    tex.close();


    ofile->Close();

    return;
}

void TriggerEffCalib__Signal__2D(const char *labelSIGNAL, const char *inputfileSIGNAL, const char *SignalLabel) {

    gErrorIgnoreLevel = kError;

    TFile *ifileSIGNAL = new TFile(inputfileSIGNAL, "READ");

    if (!ifileSIGNAL || ifileSIGNAL->IsZombie()) {
        std::cerr << "Error: Could not open input file " << inputfileSIGNAL << std::endl;
        return;
    }

    // ------------------------------------------------------------------
    // configuration
    // ------------------------------------------------------------------
    // observables (the "VS_PseudoMET" 2D maps) and, for each, the list of
    // trigger tags whose efficiency we want w.r.t. the inclusive denominator.
    struct Observable {
        std::string name;                  // denominator histo suffix
        std::vector<std::string> triggers; // numerator trigger tags
    };

    std::vector<std::string> baseTriggers = {
        "HLT_PFMET120_PFMHT120_IDTight",
        "HLT_PFHT500_PFMET100_PFMHT100_IDTight",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",
        "HLT_MET105_IsoTrk50",
        "orMETtrg"
    };

    std::vector<std::string> orMET3a4 = {
        "orMET3a4trg1", "orMET3a4trg2", "orMET3a4trg3", "orMET3a4trg4"
    };

    // PuppiMET_VS_PseudoMET : original triggers + the 4 orMET3a4trg
    std::vector<std::string> trg_PuppiMET = baseTriggers;
    trg_PuppiMET.insert(trg_PuppiMET.end(), orMET3a4.begin(), orMET3a4.end());

    // PUppiMETNoMu_VS_PseudoMET : same base triggers + the 4 orMET3a4trg
    std::vector<std::string> trg_PUppiMETNoMu = baseTriggers;
    trg_PUppiMETNoMu.insert(trg_PUppiMETNoMu.end(), orMET3a4.begin(), orMET3a4.end());

    std::vector<Observable> observables = {
        { "PUppiMET_VS_PseudoMET",     trg_PuppiMET     },
        { "PUppiMETNoMu_VS_PseudoMET", trg_PUppiMETNoMu }
    };

    const char *xtitle = "PUppi MET [GeV]";
    const char *ytitle = "Pseudo MET [GeV]";

    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);

    // ------------------------------------------------------------------
    // loop over observables / triggers
    // ------------------------------------------------------------------
    for (const auto &obs : observables) {

        const std::string &obsName = obs.name;

        // inclusive denominator
        TH2F *den = (TH2F*)ifileSIGNAL->Get(Form("%s_%s", labelSIGNAL, obsName.c_str()));

        if (!den) {
            std::cerr << "Warning: missing denominator for " << obsName
                      << ", skipping." << std::endl;
            continue;
        }

        if (obsName == "PUppiMETNoMu_VS_PseudoMET") xtitle = "PUppi MET (NoMu) [GeV]";

        for (const std::string &trg : obs.triggers) {

            // numerator histo suffix: if___<trigger>___<observable>
            std::string numSuffix = "if___" + trg + "___" + obsName;

            TH2F *num = (TH2F*)ifileSIGNAL->Get(Form("%s_%s", labelSIGNAL, numSuffix.c_str()));

            if (!num) {
                std::cerr << "Warning: missing numerator " << numSuffix
                          << ", skipping." << std::endl;
                continue;
            }

            // ---- efficiency = numerator / denominator ----
            TH2F *eff = (TH2F*)num->Clone(Form("eff_%s_%s_SIGNAL", trg.c_str(), obsName.c_str()));
            eff->SetDirectory(0);
            eff->Divide(den);

            // ---- drawing ----
            TCanvas *c = DrawCanvas(eff, Form("eff_%s_%s_SIGNAL", trg.c_str(), obsName.c_str()),
                                    xtitle, ytitle, Form("Eff. %s", trg.c_str()),
                                    "COLZ", 0, 1200, 0, 1200);
            c->cd(); latex1->Draw();

            // ---- saving ----
            c->SaveAs(Form("../outputDisplay/TriggEff/c_%s_%s_SIGNAL__%s.pdf",
                           trg.c_str(), obsName.c_str(), SignalLabel));
        }
    }

    return;
}


void TriggerEffCalib (const char *labelData, const char *labelMC,
                      const char *inputfileDATA, const char *inputfileMC,
                      const char *ofilename) {

    gErrorIgnoreLevel = kError;

    cout << "dataset: " << inputfileDATA << " and " << inputfileMC << endl;

    TFile *ifileDATA = new TFile(inputfileDATA, "READ");
    TFile *ifileMC   = new TFile(inputfileMC,   "READ");

    if (!ifileDATA || ifileDATA->IsZombie()) {
        std::cerr << "Error: Could not open input file " << inputfileDATA << std::endl;
        return;
    }
    if (!ifileMC || ifileMC->IsZombie()) {
        std::cerr << "Error: Could not open input file " << inputfileMC << std::endl;
        return;
    }

    // ------------------------------------------------------------------
    // configuration
    // ------------------------------------------------------------------
    struct Observable {
        std::string labelData;             // selection label to read in DATA
        std::string labelMC;               // selection label to read in MC
        std::string branch;                // histo suffix in the file
        std::string tag;                   // short tag for canvas / PDF names
        std::string xtitle;                // x axis title
        double cut;                        // lower threshold [GeV] for the "after-cut" efficiency (<0 = none)
        std::vector<std::string> triggers; // numerator trigger tags
    };

    std::vector<std::string> baseTriggers = {
        "HLT_PFMET120_PFMHT120_IDTight",
        "HLT_PFHT500_PFMET100_PFMHT100_IDTight",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",
        "HLT_MET105_IsoTrk50",
        "orMETtrg"
    };

    std::vector<std::string> orMET3a4 = {
        "orMET3a4trg1", "orMET3a4trg2", "orMET3a4trg3", "orMET3a4trg4"
    };

    // base triggers + the 4 orMET3a4trg
    std::vector<std::string> fullTriggers = baseTriggers;
    fullTriggers.insert(fullTriggers.end(), orMET3a4.begin(), orMET3a4.end());

    // marker sizes used in the original macro, indexed by base trigger
    std::map<std::string, double> markerSize = {
        { "HLT_PFMET120_PFMHT120_IDTight",                 0.045 },
        { "HLT_PFHT500_PFMET100_PFMHT100_IDTight",         0.04  },
        { "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",  0.032 },
        { "HLT_MET105_IsoTrk50",                           0.05  },
        { "orMETtrg",                                      0.06  }
    };

    std::string lD  = labelData;
    std::string lMC = labelMC;

    std::vector<Observable> observables = {
        // labelData, labelMC, branch,          tag,             xtitle,                   cut,   triggers
        { lD,  lMC,  "PseudoCaloMET", "PseudoMET",     "Pseudo MET [GeV]",       250.,  fullTriggers },
        { lD,  lMC,  "PUppiMET",      "PUppiMET",      "PUppi MET [GeV]",        150.,  fullTriggers },
        { lD,  lMC,  "PUppiMETNoMu",  "PUppiMETNoMu",  "PUppi MET (NoMu) [GeV]", 150.,  fullTriggers },
        { lD,  lMC,  "PFtrackPT",     "PFtrackPT",     "PF track p_{T} [GeV]",   -1.,   { "HLT_MET105_IsoTrk50" } }
    };

    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);

    // fichiers SF par bin pour orMETtrg (PseudoCaloMET et PUppiMET)
    std::ofstream sfPseudo(Form("../outputDisplay/TriggEff/SF_orMETtrg_PseudoCaloMET_%s.txt", ofilename));
    std::ofstream sfPUppi (Form("../outputDisplay/TriggEff/SF_orMETtrg_PUppiMET_%s.txt",      ofilename));

    std::ostringstream texBody;

    // binomial efficiency + error from integrated counts.
    // counts are summed over bins with low edge >= xmin (xmin < 0 -> full range,
    // under/overflow included).
    auto integratedEff = [](TH1F *num, TH1F *den, double xmin,
                            double &eff, double &err) {
        int blo, bhi;
        if (xmin < 0) {
            blo = 0;                       // include underflow
            bhi = den->GetNbinsX() + 1;    // include overflow
        } else {
            blo = den->GetXaxis()->FindBin(xmin);
            bhi = den->GetNbinsX() + 1;    // up to overflow
        }
        double N = den->Integral(blo, bhi);
        double k = num->Integral(blo, bhi);
        if (N > 0) {
            eff = k / N;
            double var = eff * (1.0 - eff) / N;   // gaussian binomial approx
            err = (var > 0) ? std::sqrt(var) : 0.0;
        } else { eff = 0.0; err = 0.0; }
    };

    // ------------------------------------------------------------------
    // loop over observables / triggers
    // ------------------------------------------------------------------
    for (const auto &obs : observables) {

        // inclusive denominators
        TH1F *den_DATA = (TH1F*)ifileDATA->Get(Form("%s_%s", obs.labelData.c_str(), obs.branch.c_str()));
        TH1F *den_MC   = (TH1F*)ifileMC->Get(Form("%s_%s",   obs.labelMC.c_str(),   obs.branch.c_str()));

        if (!den_DATA || !den_MC) {
            std::cerr << "Warning: missing denominator for " << obs.branch
                      << " (DATA=" << den_DATA << ", MC=" << den_MC << "), skipping observable." << std::endl;
            continue;
        }

        double xUp = (obs.branch == "PFtrackPT") ? 1500 : 1200;

        for (const std::string &trg : obs.triggers) {

            std::string numName_DATA = Form("%s_if___%s___%s",
                                            obs.labelData.c_str(), trg.c_str(), obs.branch.c_str());
            std::string numName_MC   = Form("%s_if___%s___%s",
                                            obs.labelMC.c_str(),   trg.c_str(), obs.branch.c_str());

            TH1F *num_DATA = (TH1F*)ifileDATA->Get(numName_DATA.c_str());
            TH1F *num_MC   = (TH1F*)ifileMC->Get(numName_MC.c_str());

            if (!num_DATA || !num_MC) {
                std::cerr << "Warning: missing numerator " << trg << " / " << obs.branch
                          << " (DATA=" << num_DATA << ", MC=" << num_MC << "), skipping." << std::endl;
                continue;
            }

            if (obs.branch == "PFtrackPT") { num_MC->Rebin(3); num_DATA->Rebin(3); den_MC->Rebin(3); den_DATA->Rebin(3); }

            // ---- efficiency = numerator / denominator ----
            TH1F *eff_DATA = (TH1F*)num_DATA->Clone(Form("eff_%s_%s_DATA", trg.c_str(), obs.tag.c_str()));
            TH1F *eff_MC   = (TH1F*)num_MC->Clone(Form("eff_%s_%s_MC",     trg.c_str(), obs.tag.c_str()));
            eff_DATA->SetDirectory(0);
            eff_MC->SetDirectory(0);
            eff_DATA->Divide(den_DATA);
            eff_MC->Divide(den_MC);

            // ---- écriture du SF par bin pour orMETtrg ----
            if (trg == "orMETtrg" &&
                (obs.branch == "PseudoCaloMET" || obs.branch == "PUppiMET")) {

                std::ofstream &out = (obs.branch == "PseudoCaloMET") ? sfPseudo : sfPUppi;

                for (int b = 1; b <= eff_DATA->GetNbinsX(); ++b) {
                    double eD  = eff_DATA->GetBinContent(b);
                    double eM  = eff_MC->GetBinContent(b);
                    double erD = eff_DATA->GetBinError(b);
                    double erM = eff_MC->GetBinError(b);

                    double sf = (eM > 0) ? eD / eM : 0.0;

                    // propagation d'erreur du quotient : (dSF/SF)^2 = (dD/D)^2 + (dM/M)^2
                    double sfErr = 0.0;
                    if (eM > 0 && eD > 0) {
                        double relD = erD / eD;
                        double relM = erM / eM;
                        sfErr = sf * std::sqrt(relD * relD + relM * relM);
                    } else if (eM > 0) {
                        sfErr = erD / eM;   // cas eD = 0
                    }

                    double sfDown = sf - sfErr;
                    double sfUp   = sf + sfErr;
                    double xup    = eff_DATA->GetXaxis()->GetBinUpEdge(b);

                    out << xup << " " << sfDown << " " << sf << " " << sfUp << "\n";
                }
            }

            // ---- integrated efficiencies for the LaTeX table ----
            // full range
            double effD = 0., errD = 0., effM = 0., errM = 0.;
            integratedEff(num_DATA, den_DATA, -1., effD, errD);
            integratedEff(num_MC,   den_MC,   -1., effM, errM);

            // after a lower cut on the observable (if defined for this obs)
            double effDc = -1., errDc = 0., effMc = -1., errMc = 0., sfc = 0., sfcErr = 0.;
            if (obs.cut >= 0) {
                integratedEff(num_DATA, den_DATA, obs.cut, effDc, errDc);
                integratedEff(num_MC,   den_MC,   obs.cut, effMc, errMc);
            }

            std::string obsTex = obs.tag; std::string trgTex = trg;
            for (auto *s : { &obsTex, &trgTex })
                for (size_t p = 0; (p = s->find('_', p)) != std::string::npos; p += 2)
                    s->replace(p, 1, "\\_");

            // after-cut cells: values if a cut is defined, dashes otherwise
            std::string cellDc, cellMc, cellSFc;
            if (obs.cut >= 0) {
                cellDc  = Form("%.4f $\\pm$ %.4f", effDc, errDc);
                cellMc  = Form("%.4f $\\pm$ %.4f", effMc, errMc);
            } else {
                cellDc = cellMc = cellSFc = "--";
            }

            texBody << "\\texttt{" << trgTex << "} & \\texttt{" << obsTex << "} & "
                    << cellDc << " & " << cellMc << " \\\\\n";

            // ---- drawing : DATA + MC superimposed ----
            double msize = markerSize.count(trg) ? markerSize[trg] : 0.045;
            std::string cName  = Form("c_%s___%s", trg.c_str(), obs.tag.c_str());
            std::string ytitle = (trg == "orMETtrg") ? "eff orMETtrg"
                                                      : Form("eff. %s", trg.c_str());

            latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Work in progress}");
            TCanvas *c = DrawCanvas(eff_DATA, eff_MC, cName.c_str(),
                                    obs.xtitle.c_str(), ytitle.c_str(),
                                    "E1", "E1 same", true, "DATA", "MC", "lep", "lep",
                                    msize, 0, xUp, 0, 1, false);
            c->cd();
            latex1->Draw();

            // ---- ratio DATA/MC ----
            std::string ratioTag = Form("%s_%s", trg.c_str(), obs.tag.c_str());

            TCanvas *cRatio = DrawWithRatio(eff_DATA, eff_MC, c,
                                            ratioTag.c_str(), "DATA/MC",
                                            obs.xtitle.c_str(), 0, xUp);
            cRatio->SaveAs(Form("../outputDisplay/TriggEff/cRatio_%s_%s.pdf",
                                ratioTag.c_str(), ofilename));

            // ---- "bis" pass: refresh canvas then redraw (fixes rendering) ----
            latex1->SetTitle("#it{Private work (CMS simulation/data)}");
            cRatio->Modified();
            cRatio->Update();
            TCanvas *cRatio_bis = DrawWithRatio(eff_DATA, eff_MC, c,
                                                ratioTag.c_str(), "DATA/MC",
                                                obs.xtitle.c_str(), 0, xUp);
            cRatio_bis->SaveAs(Form("../outputDisplay/TriggEff/cRatio_%s_%s_bis.pdf",
                                    ratioTag.c_str(), ofilename));
            cRatio_bis->SaveAs(Form("../outputDisplay/TriggEff/cRatio_%s_%s_bis.C",
                                    ratioTag.c_str(), ofilename));
        }
    }

    sfPseudo.close();
    sfPUppi.close();

    std::ofstream tex(Form("../outputDisplay/TriggEff/TriggerEff_table_%s.txt", ofilename));
    tex << "\\begin{table}[htbp]\n  \\centering\n"
        << "  \\caption{Integrated trigger efficiency, data vs.\\ MC. "
        << "The last three columns give the efficiencies and scale factor after a lower cut "
        << "on the observable (Pseudo MET $>$ 250 GeV, PUppi MET $>$ 150 GeV).}\n"
        << "  \\begin{tabular}{llcccc}\n    \\hline\n"
        << "     & \\multicolumn{2}{c}{after cut} \\\\\n"
        << "    Trigger & Observable &"
        << "$\\varepsilon_{\\mathrm{data}}$ & $\\varepsilon_{\\mathrm{MC}}$\\\\\n    \\hline\n"
        << texBody.str()
        << "    \\hline\n  \\end{tabular}\n\\end{table}\n";
    tex.close();

    return;
}

void TriggerEffCalib__2D(const char *labelData, const char *labelMC,
                         const char *inputfileDATA, const char *inputfileMC,
                         const char *ofilename) {
 
    gErrorIgnoreLevel = kError;
 
    cout << "dataset: " << inputfileDATA << " and " << inputfileMC << endl;
 
    TFile *ifileDATA = new TFile(inputfileDATA, "READ");
    TFile *ifileMC   = new TFile(inputfileMC,   "READ");
 
    if (!ifileDATA || ifileDATA->IsZombie()) {
        std::cerr << "Error: Could not open input file " << inputfileDATA << std::endl;
        return;
    }
    if (!ifileMC || ifileMC->IsZombie()) {
        std::cerr << "Error: Could not open input file " << inputfileMC << std::endl;
        return;
    }
 
    // ------------------------------------------------------------------
    // custom binning  (x = PUppi MET, y = Pseudo MET)
    // ------------------------------------------------------------------
    std::vector<double> xBins = {0, 150, 200, 300, 400, 500, 600, 700, 800, 1200};
    std::vector<double> yBins = {0, 250, 300, 400, 500, 600, 700, 800, 1200};
 
    // ------------------------------------------------------------------
    // configuration : observables + their trigger lists
    // ------------------------------------------------------------------
    struct Observable {
        std::string name;                  // denominator histo suffix
        std::vector<std::string> triggers; // numerator trigger tags
    };
 
    std::vector<std::string> baseTriggers = {
        "HLT_PFMET120_PFMHT120_IDTight",
        "HLT_PFHT500_PFMET100_PFMHT100_IDTight",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",
        "HLT_MET105_IsoTrk50",
        "orMETtrg"
    };
 
    std::vector<std::string> orMET3a4 = {
        "orMET3a4trg1", "orMET3a4trg2", "orMET3a4trg3", "orMET3a4trg4"
    };
 
    // base triggers + the 4 orMET3a4trg (same list for both observables)
    std::vector<std::string> triggers = baseTriggers;
    triggers.insert(triggers.end(), orMET3a4.begin(), orMET3a4.end());
 
    std::vector<Observable> observables = {
        { "PUppiMET_VS_PseudoMET",     triggers },
        { "PUppiMETNoMu_VS_PseudoMET", triggers }
    };
 
    const char *xtitle = "PUppi MET [GeV]";
    const char *ytitle = "Pseudo MET [GeV]";
 
    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Private work (CMS simulation/data)}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);
 
    // ------------------------------------------------------------------
    // loop over observables
    // ------------------------------------------------------------------
    for (const auto &obs : observables) {
 
        const std::string &obsName = obs.name;
 
        // adapt x-axis title for the NoMu map
        const char *xtitle_obs = (obsName == "PUppiMETNoMu_VS_PseudoMET")
                               ? "PUppi MET (NoMu) [GeV]" : xtitle;
 
        // --------------------------------------------------------------
        // inclusive denominators
        // --------------------------------------------------------------
        TH2F *den_MC_raw   = (TH2F*)ifileMC->Get(Form("%s_%s",   labelMC,   obsName.c_str()));
        TH2F *den_DATA_raw = (TH2F*)ifileDATA->Get(Form("%s_%s", labelData, obsName.c_str()));
 
        if (!den_MC_raw || !den_DATA_raw) {
            std::cerr << "Error: missing denominator for " << obsName
                      << " (MC=" << den_MC_raw << ", DATA=" << den_DATA_raw << "), skipping." << std::endl;
            continue;
        }
 
        // rebinned denominators (once per observable)
        TH2D *den_MC_rb   = RebinTH2(den_MC_raw,   xBins, yBins, Form("den_MC_%s_rebin",   obsName.c_str()));
        TH2D *den_DATA_rb = RebinTH2(den_DATA_raw, xBins, yBins, Form("den_DATA_%s_rebin", obsName.c_str()));
 
        // --------------------------------------------------------------
        // loop over triggers
        // --------------------------------------------------------------
        for (const std::string &trg : obs.triggers) {
 
            // numerator histo suffix: if___<trigger>___<observable>
            std::string numSuffix = "if___" + trg + "___" + obsName;
 
            TH2F *num_MC_raw   = (TH2F*)ifileMC->Get(Form("%s_%s",   labelMC,   numSuffix.c_str()));
            TH2F *num_DATA_raw = (TH2F*)ifileDATA->Get(Form("%s_%s", labelData, numSuffix.c_str()));
 
            if (!num_MC_raw || !num_DATA_raw) {
                std::cerr << "Warning: missing numerator " << numSuffix
                          << " (MC=" << num_MC_raw << ", DATA=" << num_DATA_raw << "), skipping." << std::endl;
                continue;
            }
 
            // ==========================================================
            // (A) FINE BINNING : efficiency / SF on the original binning
            // ==========================================================
            TH2F *eff_MC_fine   = (TH2F*)num_MC_raw->Clone(Form("eff_%s_%s_MC_fineBin",   trg.c_str(), obsName.c_str()));
            TH2F *eff_DATA_fine = (TH2F*)num_DATA_raw->Clone(Form("eff_%s_%s_DATA_fineBin", trg.c_str(), obsName.c_str()));
            eff_MC_fine->SetDirectory(0);
            eff_DATA_fine->SetDirectory(0);
            eff_MC_fine->Divide(den_MC_raw);
            eff_DATA_fine->Divide(den_DATA_raw);
 
            TH2F *sf_fine = (TH2F*)eff_DATA_fine->Clone(Form("sf_%s_%s_fineBin", trg.c_str(), obsName.c_str()));
            sf_fine->SetDirectory(0);
            sf_fine->Divide(eff_MC_fine);
 
            // ==========================================================
            // (B) REBINNED : efficiency / SF on the custom binning
            // ==========================================================
            TH2D *eff_MC_rb = RebinTH2(num_MC_raw, xBins, yBins,
                                       Form("eff_%s_%s_MC_rebin", trg.c_str(), obsName.c_str()));
            TH2D *eff_DATA_rb = RebinTH2(num_DATA_raw, xBins, yBins,
                                         Form("eff_%s_%s_DATA_rebin", trg.c_str(), obsName.c_str()));
            eff_MC_rb->SetDirectory(0);
            eff_DATA_rb->SetDirectory(0);
            eff_MC_rb->Divide(den_MC_rb);
            eff_DATA_rb->Divide(den_DATA_rb);
 
            TH2D *sf_rb = (TH2D*)eff_DATA_rb->Clone(Form("sf_%s_%s_rebin", trg.c_str(), obsName.c_str()));
            sf_rb->SetDirectory(0);
            sf_rb->Divide(eff_MC_rb);
 
            // ==========================================================
            // drawing + saving : fine binning
            // ==========================================================
            TCanvas *c_MC_fine = DrawCanvas(eff_MC_fine, Form("eff_%s_%s_MC_fineBin", trg.c_str(), obsName.c_str()),
                                            xtitle_obs, ytitle, Form("Eff. %s (MC)", trg.c_str()),
                                            "COLZ", 0, 1200, 0, 1200);
            c_MC_fine->cd(); latex1->Draw();
            TCanvas *c_DATA_fine = DrawCanvas(eff_DATA_fine, Form("eff_%s_%s_DATA_fineBin", trg.c_str(), obsName.c_str()),
                                              xtitle_obs, ytitle, Form("Eff. %s (DATA)", trg.c_str()),
                                              "COLZ", 0, 1200, 0, 1200);
            c_DATA_fine->cd(); latex1->Draw();
            TCanvas *c_SF_fine = DrawCanvas(sf_fine, Form("sf_%s_%s_fineBin", trg.c_str(), obsName.c_str()),
                                            xtitle_obs, ytitle, Form("SF (DATA/MC) %s", trg.c_str()),
                                            "COLZ", 0, 1200, 0, 1200, 0, 1);
            c_SF_fine->cd(); latex1->Draw();
 
            c_MC_fine->SaveAs(Form("../outputDisplay/TriggEff/c_%s_%s_MC_fineBin__%s.pdf",     trg.c_str(), obsName.c_str(), ofilename));
            c_DATA_fine->SaveAs(Form("../outputDisplay/TriggEff/c_%s_%s_DATA_fineBin__%s.pdf", trg.c_str(), obsName.c_str(), ofilename));
            c_SF_fine->SaveAs(Form("../outputDisplay/TriggEff/c_%s_%s_SF_fineBin__%s.pdf",     trg.c_str(), obsName.c_str(), ofilename));
 
            // ==========================================================
            // drawing + saving : rebinned
            // ==========================================================
            TCanvas *c_MC_rb = DrawCanvas(eff_MC_rb, Form("eff_%s_%s_MC_rebin", trg.c_str(), obsName.c_str()),
                                          xtitle_obs, ytitle, Form("Eff. %s (MC)", trg.c_str()),
                                          "COLZ", 0, 1200, 0, 1200);
            c_MC_rb->cd(); latex1->Draw();
            TCanvas *c_DATA_rb = DrawCanvas(eff_DATA_rb, Form("eff_%s_%s_DATA_rebin", trg.c_str(), obsName.c_str()),
                                            xtitle_obs, ytitle, Form("Eff. %s (DATA)", trg.c_str()),
                                            "COLZ", 0, 1200, 0, 1200);
            c_DATA_rb->cd(); latex1->Draw();
            TCanvas *c_SF_rb = DrawCanvas(sf_rb, Form("sf_%s_%s_rebin", trg.c_str(), obsName.c_str()),
                                          xtitle_obs, ytitle, Form("SF (DATA/MC) %s", trg.c_str()),
                                          "TEXT", 0, 1200, 0, 1200, 0, 1);
            c_SF_rb->cd(); latex1->Draw();
 
            c_MC_rb->SaveAs(Form("../outputDisplay/TriggEff/c_%s_%s_MC_rebin__%s.pdf",     trg.c_str(), obsName.c_str(), ofilename));
            c_DATA_rb->SaveAs(Form("../outputDisplay/TriggEff/c_%s_%s_DATA_rebin__%s.pdf", trg.c_str(), obsName.c_str(), ofilename));
            c_SF_rb->SaveAs(Form("../outputDisplay/TriggEff/c_%s_%s_SF_rebin__%s.pdf",     trg.c_str(), obsName.c_str(), ofilename));
        }
    }
 
    return;
}


void SignalEffVsMass(const char *ofilename = "EffVsMass") {

    gErrorIgnoreLevel = kError;

    // ------------------------------------------------------------------
    // configuration
    // ------------------------------------------------------------------
    // gluino mass points (GeV) and matching file names
    std::vector<int> masses = {1100, 1200, 1300, 1400, 1600, 1800, 2000, 2200, 2400, 2600};

    auto fileName = [&](int m) {
        return Form("../output/Gluino_V19/Gluino_Run3_MET_madgraph_%d_V19p8.root", m);
    };

    // the 3 selections -> one plot each
    std::vector<std::string> AfterHSCPsel = {
        "METanalysis_TestPseudoMETCut_Eta2p4",
        "METanalysis_TestPUppiMETCut_Eta2p4",
        "METanalysis_TestPseudoMETCut_TestPUppiMETCut_Eta2p4"
    };

    // after-trigger numerator histo name, one per selection (same index order)

    std::vector<std::string> AfterTriggerNameRaw = {
        "PostTrigger_P",
        "PostTrigger_P",
        "PostTrigger_P"
    };

    std::vector<std::string> AfterTriggerNameBis = {
        "PostTriggerbis_P",
        "PostTriggerbis_P",
        "PostTriggerbis_P"
    };

    std::vector<std::string> AfterTriggerNameCut = {
        "PostTrigger_PseudoCaloMETcut_P",
        "PostTrigger_PUppiMETcut_P",
        "PostTrigger_PseudoCaloMETcut_PUppiMETcut_P"
    };

    std::vector<std::string> AfterTriggerSFNameCut = {
        "PostTriggerWsf_PseudoCaloMETcut_P",
        "PostTriggerWsf_PUppiMETcut_P",
        "PostTriggerWsf_PseudoCaloMETcut_PUppiMETcut_P"
    };
    std::vector<std::string> AfterTriggerSFupNameCut = {
        "PostTriggerWsfup_PseudoCaloMETcut_P",
        "PostTriggerWsfup_PUppiMETcut_P",
        "PostTriggerWsfup_PseudoCaloMETcut_PUppiMETcut_P"
    };
    std::vector<std::string> AfterTriggerSFdownNameCut = {
        "PostTriggerWsfdown_PseudoCaloMETcut_P",
        "PostTriggerWsfdown_PUppiMETcut_P",
        "PostTriggerWsfdown_PseudoCaloMETcut_PUppiMETcut_P"
    };

    const std::string denName = "Nosel_P";   // common denominator histo name

    // ------------------------------------------------------------------
    // binomial integrated efficiency (full histo, under/overflow included)
    // ------------------------------------------------------------------
    auto integratedEff = [](TH1 *num, TH1 *den, double &eff, double &err) {
        double N = den->Integral(0, den->GetNbinsX() + 1);
        double k = num->Integral(0, num->GetNbinsX() + 1);
        if (N > 0) {
            eff = k / N;
            double var = eff * (1.0 - eff) / N;   // gaussian binomial approx
            err = (var > 0) ? std::sqrt(var) : 0.0;
        } else { eff = 0.0; err = 0.0; }
    };

    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);

    TLatex *tex = new TLatex(0.54, 0.91, "#bf{Cut = Pseudo MET > 250 GeV}");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);
    TLatex *tex2 = new TLatex(0.56, 0.95, "#bf{Cut = Pseudo MET > 250 GeV OR}");
    tex2->SetNDC();
    tex2->SetTextFont(42);
    tex2->SetTextSize(0.04);

    // ------------------------------------------------------------------
    // one plot per selection (iterate by index to fetch the matching label)
    // ------------------------------------------------------------------
    for (size_t isel = 0; isel < AfterHSCPsel.size(); ++isel) {

        const std::string &sel      = AfterHSCPsel[isel];

        const std::string &trigNameraw = AfterTriggerNameRaw[isel];
        const std::string &trigNamebis = AfterTriggerNameBis[isel];
        const std::string &trigName = AfterTriggerNameCut[isel];
        const std::string &trigNameSF = AfterTriggerSFNameCut[isel];
        const std::string &trigNameSFup = AfterTriggerSFupNameCut[isel];
        const std::string &trigNameSFdown = AfterTriggerSFdownNameCut[isel];

        // curves:
        //   gSel        : (sel + "_P")      / Nosel_P
        //   gSelsfup    : (sel + "_Psfup")  / Nosel_P
        //   gSelsfdown  : (sel + "_Psfdown")/ Nosel_P
        //   gTrig       : (trigName)        / Nosel_P
        TGraphErrors *gSel  = new TGraphErrors();
        TGraphErrors *gSel2  = new TGraphErrors();
        TGraphErrors *gSelsfup = new TGraphErrors();
        TGraphErrors *gSelsfdown = new TGraphErrors();
        TGraphErrors *gTrigraw = new TGraphErrors();
        TGraphErrors *gTrigbis = new TGraphErrors();
        TGraphErrors *gTrig = new TGraphErrors();
        TGraphErrors *gTrigsf = new TGraphErrors();
        TGraphErrors *gTrigsfup = new TGraphErrors();
        TGraphErrors *gTrigsfdown = new TGraphErrors();

        gSel->SetName(Form("g_%s_sel",  sel.c_str()));
        gSel2->SetName(Form("g_%s_sel2",  sel.c_str()));
        gSelsfup->SetName(Form("g_%ssfup", sel.c_str()));
        gSelsfdown->SetName(Form("g_%ssfdown", sel.c_str()));
        gTrigraw->SetName(Form("g_%s_trigraw", sel.c_str()));
        gTrigbis->SetName(Form("g_%s_trigbis", sel.c_str()));
        gTrig->SetName(Form("g_%s_trig", sel.c_str()));
        gTrigsf->SetName(Form("g_%s_trigsf", sel.c_str()));
        gTrigsfup->SetName(Form("g_%s_trigsfup", sel.c_str()));
        gTrigsfdown->SetName(Form("g_%s_trigsfdown", sel.c_str()));

        int iPoint = 0;
        for (int m : masses) {

            TFile *f = new TFile(fileName(m), "READ");
            if (!f || f->IsZombie()) {
                std::cerr << "Warning: cannot open " << fileName(m) << ", skipping mass " << m << std::endl;
                continue;
            }

            TH1 *den     = (TH1*)f->Get(denName.c_str());
            TH1 *numSel  = (TH1*)f->Get(Form("%s_P", sel.c_str()));
            TH1 *numSel2  = (TH1*)f->Get(Form("%s_P2", sel.c_str()));
            TH1 *numSelsfup   = (TH1*)f->Get(Form("%s_Psfup", sel.c_str()));
            TH1 *numSelsfdown = (TH1*)f->Get(Form("%s_Psfdown", sel.c_str()));
            TH1 *numTrigraw = (TH1*)f->Get(trigNameraw.c_str());
            TH1 *numTrigbis = (TH1*)f->Get(trigNamebis.c_str());
            TH1 *numTrig = (TH1*)f->Get(trigName.c_str());

            TH1 *numTrigSF = (TH1*)f->Get(trigNameSF.c_str());
            TH1 *numTrigSFup = (TH1*)f->Get(trigNameSFup.c_str());
            TH1 *numTrigSFdown = (TH1*)f->Get(trigNameSFdown.c_str());

            if (!den || !numSel || !numSelsfup || !numSelsfdown || !numTrig || !numTrigSF || !numTrigSFup || !numTrigSFdown) {
                std::cerr << "Warning: missing histo in " << fileName(m)
                          << " (den=" << den << ", sel=" << numSel
                          << ", selsfup=" << numSelsfup << ", selsfdown=" << numSelsfdown
                          << ", trig=" << numTrig << "), skipping mass " << m << std::endl;
                f->Close();
                continue;
            }

            double effS = 0., errS = 0., effSup = 0., errSup = 0., effSdown = 0., errSdown = 0., effT = 0., errT = 0., effTsf = 0., errTsf = 0., effTsfup = 0., errTsfup = 0., effTsfdown = 0., errTsfdown = 0.;
            double effraw = 0., errraw = 0, effbis = 0., errbis = 0;
            double effSel2 = 0., errSel2 = 0.;

            integratedEff(numSel,  den, effS, errS);
            integratedEff(numSel2,  den, effSel2, errSel2);
            integratedEff(numSelsfup,   den, effSup, errSup);
            integratedEff(numSelsfdown, den, effSdown, errSdown);
            integratedEff(numTrigraw, den, effraw, errraw);
            integratedEff(numTrigbis, den, effbis, errbis);
            integratedEff(numTrig, den, effT, errT);
            integratedEff(numTrigSF, den, effTsf, errTsf);
            integratedEff(numTrigSFup, den, effTsfup, errTsfup);
            integratedEff(numTrigSFdown, den, effTsfdown, errTsfdown);

            effS *= 0.83; errS*= 0.83;  // candidate -> event efficiency
            effSel2 *= 0.83; errSel2*= 0.83;
            effSup *= 0.83; errSup*= 0.83;
            effSdown *= 0.83; errSdown*= 0.83;
            effbis *= 0.779; errbis*= 0.779;
            effraw *= 0.779; errraw*= 0.779;
            effT *= 0.779; errT*= 0.779;
            effTsf *= 0.779; errTsf*= 0.779;
            effTsfup *= 0.779; errTsfup*= 0.779;
            effTsfdown *= 0.779; errTsfdown*= 0.779;

            gSel->SetPoint(iPoint, m, effS);
            gSel->SetPointError(iPoint, 0., errS);
            gSel2->SetPoint(iPoint, m, effSel2);
            gSel2->SetPointError(iPoint, 0., errSel2);
            gSelsfup->SetPoint(iPoint, m, effSup);
            gSelsfup->SetPointError(iPoint, 0., errSup);
            gSelsfdown->SetPoint(iPoint, m, effSdown);
            gSelsfdown->SetPointError(iPoint, 0., errSdown);
            gTrigbis->SetPoint(iPoint, m, effbis);
            gTrigbis->SetPointError(iPoint, 0., errbis);
            gTrigraw->SetPoint(iPoint, m, effraw);
            gTrigraw->SetPointError(iPoint, 0., errraw);
            gTrig->SetPoint(iPoint, m, effT);
            gTrig->SetPointError(iPoint, 0., errT);
            gTrigsf->SetPoint(iPoint, m, effTsf);
            gTrigsf->SetPointError(iPoint, 0., errTsf);
            gTrigsfup->SetPoint(iPoint, m, effTsfup);
            gTrigsfup->SetPointError(iPoint, 0., errTsfup);
            gTrigsfdown->SetPoint(iPoint, m, effTsfdown);
            gTrigsfdown->SetPointError(iPoint, 0., errTsfdown);

            ++iPoint;

            f->Close();
        }

        // ---- styling ----
        gSel->SetMarkerStyle(20);
        gSel->SetMarkerColor(kRed+1);
        gSel->SetLineColor(kRed+1);
        gSel->SetLineWidth(2);

        gSel2->SetMarkerStyle(20);
        gSel2->SetMarkerColor(kYellow+1);
        gSel2->SetLineColor(kYellow+1);
        gSel2->SetLineWidth(2);

        gSelsfup->SetMarkerStyle(22);
        gSelsfup->SetMarkerColor(kOrange+7);
        gSelsfup->SetLineColor(kOrange+7);
        gSelsfup->SetLineWidth(2);

        gSelsfdown->SetMarkerStyle(23);
        gSelsfdown->SetMarkerColor(kMagenta+1);
        gSelsfdown->SetLineColor(kMagenta+1);
        gSelsfdown->SetLineWidth(2);

        gTrigraw->SetMarkerStyle(29);
        gTrigraw->SetMarkerColor(kBlack);
        gTrigraw->SetLineColor(kBlack);
        gTrigraw->SetLineWidth(2);

        gTrigbis->SetMarkerStyle(30);
        gTrigbis->SetMarkerColor(kGray);
        gTrigbis->SetLineColor(kGray);
        gTrigbis->SetLineWidth(2);

        gTrig->SetMarkerStyle(21);
        gTrig->SetMarkerColor(kBlue+1);
        gTrig->SetLineColor(kBlue+1);
        gTrig->SetLineWidth(2);

        gTrigsf->SetMarkerStyle(24);
        gTrigsf->SetMarkerColor(kGreen+0);
        gTrigsf->SetLineColor(kGreen+0);
        gTrigsf->SetLineWidth(2);

        gTrigsfup->SetMarkerStyle(25);
        gTrigsfup->SetMarkerColor(kGreen+3);
        gTrigsfup->SetLineColor(kGreen+3);
        gTrigsfup->SetLineWidth(2);

        gTrigsfdown->SetMarkerStyle(26);
        gTrigsfdown->SetMarkerColor(kGreen-10);
        gTrigsfdown->SetLineColor(kGreen-10);
        gTrigsfdown->SetLineWidth(2);

        // ---- drawing ----
        TCanvas *c = new TCanvas(Form("c_eff_%s", sel.c_str()),
                                 Form("c_eff_%s", sel.c_str()), 800, 600);
        c->cd();
        c->SetGrid();
        c->SetLeftMargin(0.16);
        c->SetBottomMargin(0.16);

        TMultiGraph *mg = new TMultiGraph();
        mg->Add(gSel,  "PL");
        if(isel!=2) mg->Add(gSel2,  "PL");
        mg->Add(gSelsfup, "PL");
        mg->Add(gSelsfdown, "PL");
        mg->Add(gTrigraw, "PL");
        mg->Add(gTrigbis, "PL");
        mg->Add(gTrig, "PL");
        mg->Add(gTrigsf, "PL");
        mg->Add(gTrigsfup, "PL");
        mg->Add(gTrigsfdown, "PL");
        mg->Draw("A");
        mg->GetXaxis()->SetTitle("m_{#tilde{g}} [GeV]");
        mg->GetYaxis()->SetTitle("Acceptance");
        mg->SetMinimum(0.0);
        mg->SetMaximum(0.6);
        mg->GetXaxis()->SetLimits(1000, 2700);

        mg->SetTitle("");
        mg->GetYaxis()->SetTitleSize(0.06);
        mg->GetXaxis()->SetTitleSize(0.06);
        mg->GetXaxis()->SetTitleOffset(0.9);
        mg->GetYaxis()->SetTitleOffset(1);
        mg->GetXaxis()->SetLabelSize(0.05);
        mg->GetYaxis()->SetLabelSize(0.05);

        TLegend *leg = new TLegend(0.18, 0.60, 0.89, 0.89);
        leg->SetBorderSize(0);
        leg->SetNColumns(2);
        leg->AddEntry(gTrigraw, "Trigger + METfilters", "lep");
        leg->AddEntry(gTrigbis, "Trigger (no jet check) + METfilters", "lep");
        leg->AddEntry(gTrig, "Trigger + METfilters + Cut", "lep");
        if(isel!=2) leg->AddEntry(gSel2,  "HSCP pre-sel. + SF (2D)", "lep");
        if(isel==2) leg->AddEntry(gSel,  "HSCP pre-sel. + SF (2D)", "lep");
        if(isel!=2) leg->AddEntry(gTrigsf, "Trigger + METfilters + Cut + SF (1D)", "lep");
        if(isel==2) leg->AddEntry(gTrigsf, "Trigger + METfilters + Cut + SF (2D)", "lep");
        if(isel!=2) leg->AddEntry(gSel,  "HSCP pre-sel. + SF (1D)", "lep");
        if(isel==2) leg->AddEntry(gSelsfup, "HSCP pre-sel. + SF_{up} (2D)", "lep");
        if(isel!=2) leg->AddEntry(gTrigsfup, "Trigger + METfilters + Cut + SF_{up} (1D)", "lep");
        if(isel==2) leg->AddEntry(gTrigsfup, "Trigger + METfilters + Cut + SF_{up} (2D)", "lep");
        if(isel!=2) leg->AddEntry(gSelsfup, "HSCP pre-sel. + SF_{up} (1D)", "lep");
        if(isel==2) leg->AddEntry(gSelsfdown, "HSCP pre-sel. + SF_{down} (2D)", "lep");
        if(isel!=2) leg->AddEntry(gTrigsfdown, "Trigger + METfilters + Cut + SF_{down} (1D)", "lep");
        if(isel==2) leg->AddEntry(gTrigsfdown, "Trigger + METfilters + Cut + SF_{down} (2D)", "lep");
        if(isel!=2) leg->AddEntry(gSelsfdown, "HSCP pre-sel. + SF_{down} (1D)", "lep");
        leg->Draw();

        latex1->Draw();
        if (isel==1) {tex->SetTitle("#bf{Cut = PUppi MET > 150 GeV}");tex->SetX(0.56);}
        tex->Draw();
        if (isel==2) tex2->Draw();

        c->SaveAs(Form("../outputDisplay/TriggEff/c_EffVsMass_%s__%s.pdf", sel.c_str(), ofilename));

        latex1->SetTitle("#it{Private work (CMS simulation)}");
        c->Modified();
        c->Update();
        c->SaveAs(Form("../outputDisplay/TriggEff/c_EffVsMass_%s__%s_bis.pdf", sel.c_str(), ofilename));
        latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    }

    return;
}


void DisplayTriggerEff(const TString& inputFile = "../output/Gluino_V19/Gluino_Run3_MET_madgraph_2000_V19p6.root",
                       bool sig = true) {

    gStyle->SetOptStat(0);

    TFile* f = TFile::Open(inputFile);
    if (!f || f->IsZombie()) {
        std::cerr << "Error: cannot open file " << inputFile << std::endl;
        return;
    }

    const char* names[2] = {
        "TriggerEffCalib__Signal_if___orMETtrg___PseudoCaloMET",
        "TriggerEffCalib__Signal_if___orMETtrg___PUppiMET"
    };

    if (!sig) {
        names[0] = "TriggerEffCalib_if___orMETtrg___PseudoCaloMET";
        names[1] = "TriggerEffCalib_if___orMETtrg___PUppiMET";
    }

    const char* labels[2] = { "Pseudo MET", "PUppi MET" };
    const int colors[2] = { kGreen+3, kOrange+8 };

    TH1* h[2] = { nullptr, nullptr };
    for (int i = 0; i < 2; ++i) {
        h[i] = dynamic_cast<TH1*>(f->Get(names[i]));
        if (!h[i]) {
            std::cerr << "Error: histogram " << names[i] << " not found" << std::endl;
            return;
        }
        h[i]->SetDirectory(nullptr);
        h[i]->SetLineColor(colors[i]);
        h[i]->SetLineWidth(2);
        h[i]->SetMarkerColor(colors[i]);
        h[i]->SetMarkerStyle(20 + i);
    }
    f->Close();

    // ---------------- Canvas & pads ----------------
    TCanvas* c = new TCanvas("c_TriggerEff", "Trigger efficiency (orMETtrg)", 800, 600);

    TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0);
    pad1->SetLeftMargin(0.16);
    pad1->SetBottomMargin(0.16);
    pad1->Draw();

    TPad* pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.315);
    pad2->SetLeftMargin(0.16);pad2->SetBottomMargin(0.33);
    pad2->Draw();

    // ---------------- Upper pad: distributions ----------------
    pad1->cd();

    pad1->SetLogy();

    double ymax = 0;
    for (int i = 0; i < 2; ++i) ymax = std::max(ymax, h[i]->GetMaximum());
    h[0]->SetMaximum(1.2 * ymax);
    h[0]->SetMinimum(0.1);

    h[0]->SetTitle(";MET [GeV];Nb of Events");
    h[0]->GetYaxis()->SetTitleSize(0.065);
    h[0]->GetYaxis()->SetTitleOffset(0.7);
    h[0]->GetYaxis()->SetLabelSize(0.055);
    h[0]->GetXaxis()->SetLabelSize(0.055);   // labels only on lower pad
    h[0]->GetXaxis()->SetTitleSize(0.065);

    for (int i = 0; i < 2; ++i)
        h[i]->Draw(i == 0 ? "HIST" : "HIST SAME");

    TLegend* leg = new TLegend(0.60, 0.68, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    for (int i = 0; i < 2; ++i) leg->AddEntry(h[i], labels[i], "l");
    leg->Draw();

    TLatex* latex1 = new TLatex(0.16, 0.91, "#it{Private Work (CMS simulation)}");
    if (!sig) latex1->SetTitle("#it{Private Work (CMS data)}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.05);
    latex1->Draw();

    TLatex* latex2 = new TLatex(0.78, 0.91, "#scale[1.3]{#bf{m_{#tilde{g}}=2000}}");
    latex2->SetNDC();
    latex2->SetTextFont(42);
    latex2->SetTextSize(0.05);
    if (!sig) latex2->SetTitle("109 fb^{-1} (13.6 TeV)");
    if (!sig) latex2->SetX(0.71);
    latex2->Draw();

    // ---------------- Lower pad: CDFs ----------------
    TH1* hc[2] = { nullptr, nullptr };
    TH1* hCDF[2] = { nullptr, nullptr };
    for (int i = 0; i < 2; ++i) {
        hc[i] = (TH1*)h[i]->Clone(Form("%s_clone_cdf", h[i]->GetName()));
        hc[i]->SetDirectory(nullptr);
        if (hc[i]->Integral() > 0) hc[i]->Scale(1.0 / hc[i]->Integral());
        hCDF[i] = hc[i]->GetCumulative();
        hCDF[i]->SetName(Form("hCDF%d_TriggerEff", i));
        hCDF[i]->SetDirectory(nullptr);
        hCDF[i]->SetLineColor(colors[i]);
        hCDF[i]->SetMarkerColor(colors[i]);
        hCDF[i]->SetMarkerStyle(20 + i);
        hCDF[i]->SetLineWidth(2);
    }

    double Xmin = h[0]->GetXaxis()->GetXmin();
    double Xmax = h[0]->GetXaxis()->GetXmax();

    pad2->cd();
    gPad->SetTickx(0);

    hCDF[0]->SetTitle("");
    hCDF[0]->GetYaxis()->SetTitle("CDF");
    hCDF[0]->GetXaxis()->SetTitle("MET [GeV]");
    hCDF[0]->GetYaxis()->SetRangeUser(0, 1.05);
    hCDF[0]->GetYaxis()->SetNdivisions(505);
    hCDF[0]->GetYaxis()->SetTitleFont(43);
    hCDF[0]->GetXaxis()->SetTitleFont(43);
    hCDF[0]->GetYaxis()->SetLabelFont(43);
    hCDF[0]->GetXaxis()->SetLabelFont(43);
    hCDF[0]->GetYaxis()->SetTitleSize(26);
    hCDF[0]->GetXaxis()->SetTitleSize(26);
    hCDF[0]->GetYaxis()->SetLabelSize(22);
    hCDF[0]->GetXaxis()->SetLabelSize(22);
    hCDF[0]->GetYaxis()->SetTitleOffset(1.3);
    hCDF[0]->GetXaxis()->SetTitleOffset(1.0);
    hCDF[0]->GetXaxis()->SetRangeUser(Xmin, Xmax);

    hCDF[0]->Draw("P");
    hCDF[1]->Draw("P SAME");

    TLine* line = new TLine(Xmin, 0.5, Xmax, 0.5);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->Draw("same");

    c->cd();
    c->Update();
    c->SaveAs(Form("../outputDisplay/TriggEff/METdist%s.pdf", (sig) ? "_SIGNAL" : "_DATA"));
}


void FpixelInSignalAndData() {

    TFile *ifileData   = new TFile("../output/JetMET2024_V12/JetMET2024_V12p24.root", "READ");
    TFile *ifileSignal = new TFile("../output/Gluino_V19/Gluino_Run3_MET_madgraph_2000_V19p6.root", "READ");

    TH1D *FpixSignal = (TH1D*)ifileSignal->Get("METanalysis_TestPUppiMETCut_Eta2p4_Fpix");
    TH2D *FpixData_A3fp9 = (TH2D*)ifileData->Get("fpix_eta_regionA_3fp9_METanalysis_Eta2p4");
    TH2D *FpixData_A9fp10 = (TH2D*)ifileData->Get("fpix_eta_regionA_9fp10_METanalysis_Eta2p4");
    TH2D *FpixData_C3fp9 = (TH2D*)ifileData->Get("fpix_eta_regionC_3fp9_METanalysis_Eta2p4");
    TH2D *FpixData_C9fp10 = (TH2D*)ifileData->Get("fpix_eta_regionC_9fp10_METanalysis_Eta2p4");

    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);


    FpixSignal->SetBinContent(FpixSignal->GetNbinsX()-1, FpixSignal->GetBinContent(FpixSignal->GetNbinsX()-1) + FpixSignal->GetBinContent(FpixSignal->GetNbinsX()));
    FpixSignal->SetBinContent(FpixSignal->GetNbinsX(), 0);
    FpixSignal->SetBinContent(FpixSignal->GetNbinsX()-2, FpixSignal->GetBinContent(FpixSignal->GetNbinsX()-2) + FpixSignal->GetBinContent(FpixSignal->GetNbinsX()-1));
    FpixSignal->SetBinContent(FpixSignal->GetNbinsX()-1, 0);

    //
    FpixData_A9fp10->Add(FpixData_A3fp9);

    TH1D *FpixData  = FpixData_A9fp10->ProjectionX("FpixData");

    TCanvas *cSignal = DrawCanvas(FpixSignal, "", "F_{pixel}", "Entries", "E1", 0.06, 0, 0.99, 633, 1, -1, true);


    TCanvas *cSignalCDF = DrawWithCDF(FpixSignal, cSignal, "FpixelInSignal_CDF", "F_{pixel}", "Signal", 0, 0.99, 633);
    TPad* p1 = (TPad*)cSignalCDF->GetPrimitive("pad1");
    p1->cd();
    latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    latex1->Draw();
    p1->Modified(); p1->Update();
    cSignalCDF->SaveAs("../outputDisplay/FpixelInSignal.pdf");
    latex1->SetTitle("#it{Private work (CMS simulation)}");
    p1->Modified(); p1->Update();
    cSignalCDF->SaveAs("../outputDisplay/FpixelInSignal_bis.pdf");


    TCanvas *cData = DrawCanvas(FpixData, "", "F_{pixel}", "Entries", "E1", 0.06, 0, 1, 602, 1, 3000, true);

    TCanvas *cDataCDF = DrawWithCDF(FpixData, cData, "FpixelInData_CDF", "F_{pixel}", "Data", 0, 1, 602);
    TPad* p2 = (TPad*)cDataCDF->GetPrimitive("pad1");
    p2->cd();
    latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    latex1->Draw();
    p2->Modified(); p2->Update();
    cDataCDF->SaveAs("../outputDisplay/FpixelInData.pdf");
    latex1->SetTitle("#it{Private work (CMS simulation)}");
    p2->Modified(); p2->Update();
    cDataCDF->SaveAs("../outputDisplay/FpixelInData_bis.pdf");


    return;
}


void PairTypeStages_SingleMass(int mass = 2000,
                               const char *ofilename = "PairTypeStages") {

    gErrorIgnoreLevel = kError;
    gStyle->SetOptStat(0);

    auto fileName = [&](int m) {
        return Form("../output/Gluino_V19/Gluino_Run3_MET_madgraph_%d_V19p6.root", m);
    };

    // convention : valeur pairType=v -> bin ROOT v+1  => 1,2,3 dans les bins 2,3,4
    std::vector<int>         pairBins = {2, 3, 4};
    std::vector<std::string> catLabels = {"charged-charged", "neutral-charged", "neutral-neutral"};

    // les 5 stades : nom d'histo + label légende + style
    std::vector<std::string> stageHistos = {
        "Nosel_Gen__PairType",
        "Nosel_Gen__PairType__if_ORtrigger",
        "Nosel_GenHSCPmatching__PairType__0HSCP_if_ORtrigger",
        "Nosel_GenHSCPmatching__PairType__1HSCP_if_ORtrigger",
        "Nosel_GenHSCPmatching__PairType__2HSCP_if_ORtrigger"
    };
    std::vector<std::string> stageLabels = {
        "No trigger",
        "with OR trigger",
        "with OR trigger + 0 HSCP matched",
        "with OR trigger + 1 HSCP matched",
        "with OR trigger + 2 HSCP matched"
    };
    std::vector<int> stageColors  = {kBlack, kRed+1, kAzure+1, kGreen+2, kMagenta+1};
    std::vector<int> stageMarkers = {20, 21, 22, 23, 29};

    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);

    TFile *f = new TFile(fileName(mass), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: cannot open " << fileName(mass) << std::endl;
        return;
    }

    // dénominateur commun = intégrale des 3 catégories de Nosel_Gen__PairType
    TH1 *hDen = (TH1*)f->Get("Nosel_Gen__PairType");
    if (!hDen) {
        std::cerr << "Error: missing Nosel_Gen__PairType in " << fileName(mass) << std::endl;
        f->Close();
        return;
    }
    double tot = 0.;
    for (int b : pairBins) tot += hDen->GetBinContent(b);
    if (tot <= 0) {
        std::cerr << "Error: null integral of Nosel_Gen__PairType" << std::endl;
        f->Close();
        return;
    }

    TCanvas *c = new TCanvas("c_PairTypeStages", "c_PairTypeStages", 800, 600);
    c->cd();
    c->SetGrid();
    c->SetLeftMargin(0.16);
    c->SetBottomMargin(0.16);

    TLegend *leg = new TLegend(0.45, 0.63, 0.88, 0.88);
    leg->SetBorderSize(0);

    std::vector<TH1F*> hStages;
    double ymax = 0.;

    for (size_t is = 0; is < stageHistos.size(); ++is) {

        TH1 *hSrc = (TH1*)f->Get(stageHistos[is].c_str());
        if (!hSrc) {
            std::cerr << "Warning: missing " << stageHistos[is]
                      << " in " << fileName(mass) << ", skipping stage." << std::endl;
            continue;
        }

        // histo condensé à 3 bins étiquetés
        TH1F *h = new TH1F(Form("hStage%zu_m%d", is, mass),
                           "", (int)pairBins.size(), 0.5, pairBins.size() + 0.5);
        h->SetDirectory(0);
        for (size_t ic = 0; ic < pairBins.size(); ++ic) {
            h->GetXaxis()->SetBinLabel(ic + 1, catLabels[ic].c_str());
            double num  = hSrc->GetBinContent(pairBins[ic]);
            double frac = num / tot;                       // normalisé par l'intégrale de Nosel_Gen__PairType
            double err  = std::sqrt(num) / tot;            // erreur poissonienne sur le numérateur
            h->SetBinContent(ic + 1, frac);
            h->SetBinError(ic + 1, err);
            ymax = std::max(ymax, frac + err);
        }

        h->SetMarkerStyle(stageMarkers[is]);
        h->SetMarkerColor(stageColors[is]);
        h->SetMarkerSize(1.5);
        h->SetLineColor(stageColors[is]);
        h->SetLineWidth(2);

        hStages.push_back(h);
        leg->AddEntry(h, stageLabels[is].c_str(), "lep");
    }
    f->Close();

    if (hStages.empty()) { delete c; return; }

    hStages[0]->SetTitle("");
    hStages[0]->GetXaxis()->SetTitle("HSCP pair type");
    hStages[0]->GetYaxis()->SetTitle("Fraction");
    hStages[0]->GetYaxis()->SetTitleSize(0.06);
    hStages[0]->GetXaxis()->SetTitleSize(0.06);
    hStages[0]->GetYaxis()->SetTitleOffset(1);
    hStages[0]->GetXaxis()->SetTitleOffset(1);
    hStages[0]->GetYaxis()->SetLabelSize(0.05);
    hStages[0]->GetXaxis()->SetLabelSize(0.06);
    hStages[0]->SetMinimum(0.0);
    hStages[0]->SetMaximum(1);

    for (size_t j = 0; j < hStages.size(); ++j)
        hStages[j]->Draw(j == 0 ? "E1" : "E1 same");
    leg->Draw();
    latex1->Draw();

    // annotation de la masse
    TLatex *mtext = new TLatex(0.18, 0.85, Form("m_{#tilde{g}} = %d GeV", mass));
    mtext->SetNDC(); mtext->SetTextFont(42); mtext->SetTextSize(0.05);
    mtext->Draw();

    c->Modified(); c->Update();
    c->SaveAs(Form("../outputDisplay/TriggEff/c_PairTypeStages_m%d__%s.pdf", mass, ofilename));

    latex1->SetTitle("#it{Private work (CMS simulation)}");
    c->Modified(); c->Update();
    c->SaveAs(Form("../outputDisplay/TriggEff/c_PairTypeStages_m%d__%s_bis.pdf", mass, ofilename));

    return;
}

void PFType_ProportionAndTrigEff(const char *ofilename = "PFTypeProp") {

    gErrorIgnoreLevel = kError;
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);

    auto fileName = [&](int m) {
        return Form("../output/Gluino_V19/Gluino_Run3_MET_madgraph_%d_V19p6.root", m);
    };

    std::vector<int> masses = {1100, 1200, 1300, 1400, 1600, 1800, 2000, 2200, 2400, 2600};

    // catégories PF : PDG id -> label
    std::vector<int>         pfPdg    = {11, 13, 211};
    std::vector<std::string> pfLabels = {"electrons", "muons", "charged pions"};

    const std::string hPFIncl = "GenHSCPmatching__PFType";
    const std::string hPFTrig = "GenHSCPmatching__PFType__if_ORtrigger";

    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);

    TCanvas *c = new TCanvas("c_PFTypeProp", "c_PFTypeProp", 800, 600);
    c->cd();
    c->SetGrid();
    c->SetLeftMargin(0.16);
    c->SetBottomMargin(0.16);

    TLegend *leg = new TLegend(0.60, 0.45, 0.88, 0.88);
    leg->SetBorderSize(0);

    std::vector<TH1F*> hProp;   // trait plein (proportion)
    std::vector<TH1F*> hEff;    // pointillé (eff OR trigger)
    std::vector<int>   hMassVal;
    double ymax = 0.;

    for (int m : masses) {

        TFile *f = new TFile(fileName(m), "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "Warning: cannot open " << fileName(m) << ", skipping mass " << m << std::endl;
            continue;
        }

        TH1 *hIncl = (TH1*)f->Get(hPFIncl.c_str());
        TH1 *hTrig = (TH1*)f->Get(hPFTrig.c_str());
        if (!hIncl || !hTrig) {
            std::cerr << "Warning: missing PFType histo in " << fileName(m)
                      << " (incl=" << hIncl << ", trig=" << hTrig << "), skipping mass " << m << std::endl;
            f->Close();
            continue;
        }

        // intégrale de reference = somme des 3 categories de GenHSCPmatching__PFType
        double totIncl = 0.;
        for (int pdg : pfPdg) totIncl += hIncl->GetBinContent(hIncl->GetXaxis()->FindBin(pdg));
        if (totIncl <= 0) { f->Close(); continue; }

        // --- h1 : proportion des PF ---
        TH1F *hp = new TH1F(Form("hPFprop_m%d", m), "", (int)pfPdg.size(), 0.5, pfPdg.size() + 0.5);
        hp->SetDirectory(0);
        // --- h2 : eff OR trigger (meme denominateur) ---
        TH1F *he = new TH1F(Form("hPFeff_m%d", m), "", (int)pfPdg.size(), 0.5, pfPdg.size() + 0.5);
        he->SetDirectory(0);

        for (size_t ic = 0; ic < pfPdg.size(); ++ic) {
            hp->GetXaxis()->SetBinLabel(ic + 1, pfLabels[ic].c_str());
            he->GetXaxis()->SetBinLabel(ic + 1, pfLabels[ic].c_str());

            double nIncl = hIncl->GetBinContent(hIncl->GetXaxis()->FindBin(pfPdg[ic]));
            double nTrig = hTrig->GetBinContent(hTrig->GetXaxis()->FindBin(pfPdg[ic]));

            double prop = nIncl / totIncl;
            double effT = nTrig / totIncl;

            hp->SetBinContent(ic + 1, prop);
            hp->SetBinError(ic + 1, std::sqrt(nIncl) / totIncl);
            he->SetBinContent(ic + 1, effT);
            he->SetBinError(ic + 1, std::sqrt(nTrig) / totIncl);

            ymax = std::max(ymax, prop + std::sqrt(nIncl) / totIncl);
        }
        f->Close();

        hProp.push_back(hp);
        hEff.push_back(he);
        hMassVal.push_back(m);
    }

    if (hProp.empty()) { delete c; return; }

    int nCurves = (int)hProp.size();

    for (int j = 0; j < nCurves; ++j) {
        int ci = (nCurves > 1)
               ? TColor::GetColorPalette((int)((j / double(nCurves - 1)) * (gStyle->GetNumberOfColors() - 1)))
               : TColor::GetColorPalette(gStyle->GetNumberOfColors() / 2);

        // trait plein : proportion
        hProp[j]->SetMarkerStyle(20);
        hProp[j]->SetMarkerColor(ci);
        hProp[j]->SetLineColor(ci);
        hProp[j]->SetLineWidth(2);
        hProp[j]->SetLineStyle(1);

        // pointillé : eff OR trigger, même couleur
        hEff[j]->SetMarkerStyle(24);
        hEff[j]->SetMarkerColor(ci);
        hEff[j]->SetLineColor(ci);
        hEff[j]->SetLineWidth(2);
        hEff[j]->SetLineStyle(2);

        leg->AddEntry(hProp[j], Form("m_{#tilde{g}} = %d GeV", hMassVal[j]), "lep");
    }

    hProp[0]->SetTitle("");
    hProp[0]->GetXaxis()->SetTitle("PF type");
    hProp[0]->GetYaxis()->SetTitle("Fraction");
    hProp[0]->GetYaxis()->SetTitleSize(0.06);
    hProp[0]->GetXaxis()->SetTitleSize(0.06);
    hProp[0]->GetYaxis()->SetTitleOffset(1);
    hProp[0]->GetXaxis()->SetTitleOffset(1);
    hProp[0]->GetYaxis()->SetLabelSize(0.05);
    hProp[0]->GetXaxis()->SetLabelSize(0.06);
    hProp[0]->SetMinimum(0.0);
    hProp[0]->SetMaximum(1.35 * ymax);
    hProp[0]->GetXaxis()->SetRangeUser(2,3);

    // dessin : d'abord tous les traits pleins, puis tous les pointillés
    hProp[0]->Draw("E1");
    for (int j = 1; j < nCurves; ++j) hProp[j]->Draw("E1 same");
    for (int j = 0; j < nCurves; ++j) hEff[j]->Draw("E1 same");

    leg->Draw();

    // légende de style (plein vs pointillé)
    TLegend *legStyle = new TLegend(0.18, 0.75, 0.45, 0.88);
    legStyle->SetBorderSize(0);
    TH1F *hSolid = new TH1F("hSolidLeg", "", 1, 0, 1); hSolid->SetLineColor(kBlack); hSolid->SetLineWidth(2); hSolid->SetLineStyle(1); hSolid->SetMarkerStyle(20);
    TH1F *hDash  = new TH1F("hDashLeg",  "", 1, 0, 1); hDash->SetLineColor(kBlack);  hDash->SetLineWidth(2);  hDash->SetLineStyle(2);  hDash->SetMarkerStyle(24);
    legStyle->AddEntry(hSolid, "No trigger", "lp");
    legStyle->AddEntry(hDash,  "with OR trigger", "lp");
    legStyle->Draw();

    latex1->Draw();

    c->Modified(); c->Update();
    c->SaveAs(Form("../outputDisplay/TriggEff/c_PFTypeProp__%s.pdf", ofilename));

    latex1->SetTitle("#it{Private work (CMS simulation)}");
    c->Modified(); c->Update();
    c->SaveAs(Form("../outputDisplay/TriggEff/c_PFTypeProp__%s_bis.pdf", ofilename));

    return;
}

void TriggerEfficiency_ByMass(const char *ofilename = "TriggerEff") {

    gErrorIgnoreLevel = kError;
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);

    auto fileName = [&](int m) {
        return Form("../output/Gluino_V19/Gluino_Run3_MET_madgraph_%d_V19p10.root", m);
    };

    std::vector<int> masses = {2000, 2400, 2600};

    // trigger : label court -> nom d'histo "passing"
    std::vector<std::string> trigLabels = {
        "PFMET120_PFMHT120", "PFHT500_PFMET100_PFMHT100",
        "PFMETNoMu120_PFMHTNoMu120", "MET105_IsoTrk50", "OR trigger", "Mu50"
    };
    std::vector<std::string> trigHistos = {
        "Nosel_Gen__PairType__if_HLT_PFMET120_PFMHT120_IDTight",
        "Nosel_Gen__PairType__if_HLT_PFHT500_PFMET100_PFMHT100_IDTight",
        "Nosel_Gen__PairType__if_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",
        "Nosel_Gen__PairType__if_HLT_MET105_IsoTrk50",
        "Nosel_Gen__PairType__if_ORtrigger",
        "Nosel_Gen__PairType__if_HLT_Mu50"
    };

    const std::string hDenomName = "Nosel_Gen__PairType";
    std::vector<int> pairTypeVals = {1, 2, 3};

    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);

    TCanvas *c = new TCanvas("c_TriggerEfficiency", "c_TriggerEfficiency", 800, 600);
    c->cd();
    c->SetGrid();
    c->SetLeftMargin(0.16);
    c->SetBottomMargin(0.19);

    TLegend *leg = new TLegend(0.18, 0.78, 0.89, 0.89);
    leg->SetBorderSize(0);
    leg->SetNColumns(3);

    std::vector<TH1F*> hEff;
    std::vector<int>   hMassVal;
    double ymax = 0.;

    for (int m : masses) {

        TFile *f = new TFile(fileName(m), "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "Warning: cannot open " << fileName(m) << ", skipping mass " << m << std::endl;
            continue;
        }

        TH1 *hDenom = (TH1*)f->Get(hDenomName.c_str());
        if (!hDenom) {
            std::cerr << "Warning: missing " << hDenomName << " in " << fileName(m) << ", skipping mass " << m << std::endl;
            f->Close();
            continue;
        }

        double totDenom = 0.;
        for (int pt : pairTypeVals) totDenom += hDenom->GetBinContent(hDenom->GetXaxis()->FindBin(pt));
        if (totDenom <= 0) { f->Close(); continue; }

        std::vector<TH1*> hNum(trigHistos.size(), nullptr);
        bool allFound = true;
        for (size_t it = 0; it < trigHistos.size(); ++it) {
            hNum[it] = (TH1*)f->Get(trigHistos[it].c_str());
            if (!hNum[it]) {
                std::cerr << "Warning: missing " << trigHistos[it] << " in " << fileName(m) << std::endl;
                allFound = false;
            }
        }
        if (!allFound) { f->Close(); continue; }

        TH1F *he = new TH1F(Form("hTrigEff_m%d", m), "", (int)trigLabels.size(), 0.5, trigLabels.size() + 0.5);
        he->SetDirectory(0);

        for (size_t it = 0; it < trigLabels.size(); ++it) {
            he->GetXaxis()->SetBinLabel(it + 1, trigLabels[it].c_str());

            double nNum = 0.;
            for (int pt : pairTypeVals) nNum += hNum[it]->GetBinContent(hNum[it]->GetXaxis()->FindBin(pt));

            double eff = nNum / totDenom;
            he->SetBinContent(it + 1, eff);
            he->SetBinError(it + 1, std::sqrt(nNum) / totDenom);

            ymax = std::max(ymax, eff + std::sqrt(nNum) / totDenom);
        }
        f->Close();

        hEff.push_back(he);
        hMassVal.push_back(m);
    }

    if (hEff.empty()) { delete c; return; }

    int nCurves = (int)hEff.size();

    for (int j = 0; j < nCurves; ++j) {
        // int ci = (nCurves > 1)
        //        ? TColor::GetColorPalette((int)((j / double(nCurves - 1)) * (gStyle->GetNumberOfColors() - 1)))
        //        : TColor::GetColorPalette(gStyle->GetNumberOfColors() / 2);

        int ci = 1; 
        if (j==0) ci = kOrange+8;
        if (j==1) ci = kViolet+1;
        if (j==2) ci = kGreen-3;

        hEff[j]->SetMarkerStyle(20+j);
        hEff[j]->SetMarkerColor(ci);
        hEff[j]->SetLineColor(ci);
        hEff[j]->SetLineWidth(2);
        hEff[j]->SetLineStyle(1);

        leg->AddEntry(hEff[j], Form("m_{#tilde{g}} = %d GeV", hMassVal[j]), "lep");
    }

    hEff[0]->SetTitle("");
    hEff[0]->GetXaxis()->SetTitle("");
    hEff[0]->GetYaxis()->SetTitle("Efficiency");
    hEff[0]->GetYaxis()->SetTitleSize(0.06);
    hEff[0]->GetXaxis()->SetTitleSize(0.06);
    hEff[0]->GetYaxis()->SetTitleOffset(1);
    hEff[0]->GetXaxis()->SetTitleOffset(1);
    hEff[0]->GetYaxis()->SetLabelSize(0.05);
    hEff[0]->GetXaxis()->SetLabelSize(0.045);
    hEff[0]->SetMinimum(0.0);
    hEff[0]->SetMaximum(0.4);

    hEff[0]->Draw("E1");
    for (int j = 1; j < nCurves; ++j) hEff[j]->Draw("E1 same");

    leg->Draw();
    latex1->Draw();

    c->Modified(); c->Update();
    c->SaveAs(Form("../outputDisplay/TriggEff/c_TriggerEfficiency__%s.pdf", ofilename));

    latex1->SetTitle("#it{Private work (CMS simulation)}");
    c->Modified(); c->Update();
    c->SaveAs(Form("../outputDisplay/TriggEff/c_TriggerEfficiency__%s_bis.pdf", ofilename));

    return;
}

void PairTypeStages_SingleMass__Notrigger(int mass = 2000,
                               const char *ofilename = "PairTypeStages") {

    gErrorIgnoreLevel = kError;
    gStyle->SetOptStat(0);

    auto fileName = [&](int m) {
        return Form("../output/Gluino_V19/Gluino_Run3_MET_madgraph_%d_V19p10.root", m);
    };

    // convention : valeur pairType=v -> bin ROOT v+1  => 1,2,3 dans les bins 2,3,4
    std::vector<int>         pairBins = {4,3,2};
    std::vector<std::string> catLabels = {"neutral-neutral", "neutral-charged", "charged-charged"};

    // les 5 stades : nom d'histo + label légende + style
    std::vector<std::string> stageHistos = {
        "Nosel_Gen__PairType",
        "Nosel_GenHSCPmatching__PairType__0HSCP",
        "Nosel_GenHSCPmatching__PairType__1HSCP",
        "Nosel_GenHSCPmatching__PairType__2HSCP"
    };
    std::vector<std::string> stageLabels = {
        "No trigger",
        "with 0 HSCP matched",
        "with 1 HSCP matched",
        "with 2 HSCP matched"
    };
    std::vector<int> stageColors  = {kBlack, kAzure+1, kGreen+2, kMagenta+1};
    std::vector<int> stageMarkers = {20, 22, 23, 29};

    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);

    TFile *f = new TFile(fileName(mass), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: cannot open " << fileName(mass) << std::endl;
        return;
    }

    // dénominateur commun = intégrale des 3 catégories de Nosel_Gen__PairType
    TH1 *hDen = (TH1*)f->Get("Nosel_Gen__PairType");
    if (!hDen) {
        std::cerr << "Error: missing Nosel_Gen__PairType in " << fileName(mass) << std::endl;
        f->Close();
        return;
    }
    double tot = 0.;
    for (int b : pairBins) tot += hDen->GetBinContent(b);
    if (tot <= 0) {
        std::cerr << "Error: null integral of Nosel_Gen__PairType" << std::endl;
        f->Close();
        return;
    }

    TCanvas *c = new TCanvas("c_PairTypeStages", "c_PairTypeStages", 800, 600);
    c->cd();
    c->SetGrid();
    c->SetLeftMargin(0.16);
    c->SetBottomMargin(0.16);

    TLegend *leg = new TLegend(0.18, 0.79, 0.89, 0.89);
    leg->SetBorderSize(0);
    leg->SetNColumns(2);

    std::vector<TH1F*> hStages;
    double ymax = 0.;

    for (size_t is = 0; is < stageHistos.size(); ++is) {

        TH1 *hSrc = (TH1*)f->Get(stageHistos[is].c_str());
        if (!hSrc) {
            std::cerr << "Warning: missing " << stageHistos[is]
                      << " in " << fileName(mass) << ", skipping stage." << std::endl;
            continue;
        }

        // histo condensé à 3 bins étiquetés
        TH1F *h = new TH1F(Form("hStage%zu_m%d", is, mass),
                           "", (int)pairBins.size(), 0.5, pairBins.size() + 0.5);
        h->SetDirectory(0);
        for (size_t ic = 0; ic < pairBins.size(); ++ic) {
            h->GetXaxis()->SetBinLabel(ic + 1, catLabels[ic].c_str());
            double num  = hSrc->GetBinContent(pairBins[ic]);
            double frac = num / tot;                       // normalisé par l'intégrale de Nosel_Gen__PairType
            double err  = std::sqrt(num) / tot;            // erreur poissonienne sur le numérateur
            h->SetBinContent(ic + 1, frac);
            h->SetBinError(ic + 1, err);
            ymax = std::max(ymax, frac + err);
        }

        h->SetMarkerStyle(stageMarkers[is]);
        h->SetMarkerColor(stageColors[is]);
        h->SetMarkerSize(1.5);
        h->SetLineColor(stageColors[is]);
        h->SetLineWidth(2);

        hStages.push_back(h);
        leg->AddEntry(h, stageLabels[is].c_str(), "lep");
    }
    f->Close();

    if (hStages.empty()) { delete c; return; }

    hStages[0]->SetTitle("");
    hStages[0]->GetXaxis()->SetTitle("HSCP pair type");
    hStages[0]->GetYaxis()->SetTitle("Fraction");
    hStages[0]->GetYaxis()->SetTitleSize(0.06);
    hStages[0]->GetXaxis()->SetTitleSize(0.06);
    hStages[0]->GetYaxis()->SetTitleOffset(1);
    hStages[0]->GetXaxis()->SetTitleOffset(1);
    hStages[0]->GetYaxis()->SetLabelSize(0.05);
    hStages[0]->GetXaxis()->SetLabelSize(0.06);
    hStages[0]->SetMinimum(0.0);
    hStages[0]->SetMaximum(0.6);

    for (size_t j = 0; j < hStages.size(); ++j)
        hStages[j]->Draw(j == 0 ? "E1" : "E1 same");
    leg->Draw();
    latex1->Draw();

    // annotation de la masse
    TLatex *mtext = new TLatex(0.68, 0.91, "#scale[1.3]{#bf{m_{#tilde{g}}=2000 GeV}}");
    mtext->SetNDC(); mtext->SetTextFont(42); mtext->SetTextSize(0.04);
    mtext->Draw();

    c->Modified(); c->Update();
    c->SaveAs(Form("../outputDisplay/TriggEff/PairTypeStages_SingleMass__Notrigger_%d__%s.pdf", mass, ofilename));

    latex1->SetTitle("#it{Private work (CMS simulation)}");
    c->Modified(); c->Update();
    c->SaveAs(Form("../outputDisplay/TriggEff/PairTypeStages_SingleMass__Notrigger_%d__%s_bis.pdf", mass, ofilename));

    return;
}

void PairTypeStages_SingleMass__Trigger(int mass = 2000,
                               const char *ofilename = "PairTypeStages") {

    gErrorIgnoreLevel = kError;
    gStyle->SetOptStat(0);

    auto fileName = [&](int m) {
        return Form("../output/Gluino_V19/Gluino_Run3_MET_madgraph_%d_V19p10.root", m);
    };

    // convention : valeur pairType=v -> bin ROOT v+1  => 1,2,3 dans les bins 2,3,4
    std::vector<int>         pairBins = {4,3,2};
    std::vector<std::string> catLabels = {"neutral-neutral", "neutral-charged", "charged-charged"};

    // les 5 stades : nom d'histo + label légende + style
    std::vector<std::string> stageHistos = {
        "Nosel_Gen__PairType",
        "Nosel_Gen__PairType__if_ORtrigger",
        "Nosel_Gen__PairType__if_HLT_Mu50"
    };
    std::vector<std::string> stageLabels = {
        "No trigger",
        "OR MET trigger",
        "Muon trigger"
    };
    std::vector<int> stageColors  = {kBlack, kRed+1, kYellow+2};
    std::vector<int> stageMarkers = {20, 21, 34};

    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);

    TFile *f = new TFile(fileName(mass), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: cannot open " << fileName(mass) << std::endl;
        return;
    }

    // dénominateur commun = intégrale des 3 catégories de Nosel_Gen__PairType
    TH1 *hDen = (TH1*)f->Get("Nosel_Gen__PairType");
    if (!hDen) {
        std::cerr << "Error: missing Nosel_Gen__PairType in " << fileName(mass) << std::endl;
        f->Close();
        return;
    }
    double tot = 0.;
    for (int b : pairBins) tot += hDen->GetBinContent(b);
    if (tot <= 0) {
        std::cerr << "Error: null integral of Nosel_Gen__PairType" << std::endl;
        f->Close();
        return;
    }

    TCanvas *c = new TCanvas("c_PairTypeStages", "c_PairTypeStages", 800, 600);
    c->cd();
    c->SetGrid();
    c->SetLeftMargin(0.16);
    c->SetBottomMargin(0.16);

    TLegend *leg = new TLegend(0.18, 0.79, 0.89, 0.89);
    leg->SetBorderSize(0);
    leg->SetNColumns(3);

    std::vector<TH1F*> hStages;
    double ymax = 0.;

    for (size_t is = 0; is < stageHistos.size(); ++is) {

        TH1 *hSrc = (TH1*)f->Get(stageHistos[is].c_str());
        if (!hSrc) {
            std::cerr << "Warning: missing " << stageHistos[is]
                      << " in " << fileName(mass) << ", skipping stage." << std::endl;
            continue;
        }

        // histo condensé à 3 bins étiquetés
        TH1F *h = new TH1F(Form("hStage%zu_m%d", is, mass),
                           "", (int)pairBins.size(), 0.5, pairBins.size() + 0.5);
        h->SetDirectory(0);
        for (size_t ic = 0; ic < pairBins.size(); ++ic) {
            h->GetXaxis()->SetBinLabel(ic + 1, catLabels[ic].c_str());
            double num  = hSrc->GetBinContent(pairBins[ic]);
            double frac = num / tot;                       // normalisé par l'intégrale de Nosel_Gen__PairType
            double err  = std::sqrt(num) / tot;            // erreur poissonienne sur le numérateur
            h->SetBinContent(ic + 1, frac);
            h->SetBinError(ic + 1, err);
            ymax = std::max(ymax, frac + err);
        }

        h->SetMarkerStyle(stageMarkers[is]);
        h->SetMarkerColor(stageColors[is]);
        h->SetMarkerSize(1.5);
        h->SetLineColor(stageColors[is]);
        h->SetLineWidth(2);

        hStages.push_back(h);
        leg->AddEntry(h, stageLabels[is].c_str(), "lep");
    }
    f->Close();

    if (hStages.empty()) { delete c; return; }

    hStages[0]->SetTitle("");
    hStages[0]->GetXaxis()->SetTitle("HSCP pair type");
    hStages[0]->GetYaxis()->SetTitle("Fraction");
    hStages[0]->GetYaxis()->SetTitleSize(0.06);
    hStages[0]->GetXaxis()->SetTitleSize(0.06);
    hStages[0]->GetYaxis()->SetTitleOffset(1);
    hStages[0]->GetXaxis()->SetTitleOffset(1);
    hStages[0]->GetYaxis()->SetLabelSize(0.05);
    hStages[0]->GetXaxis()->SetLabelSize(0.06);
    hStages[0]->SetMinimum(0.0);
    hStages[0]->SetMaximum(0.6);

    for (size_t j = 0; j < hStages.size(); ++j)
        hStages[j]->Draw(j == 0 ? "E1" : "E1 same");
    leg->Draw();
    latex1->Draw();

    // annotation de la masse
    TLatex *mtext = new TLatex(0.68, 0.91, "#scale[1.3]{#bf{m_{#tilde{g}}=2000 GeV}}");
    mtext->SetNDC(); mtext->SetTextFont(42); mtext->SetTextSize(0.04);
    mtext->Draw();

    c->Modified(); c->Update();
    c->SaveAs(Form("../outputDisplay/TriggEff/PairTypeStages_SingleMass__Trigger_%d__%s.pdf", mass, ofilename));

    latex1->SetTitle("#it{Private work (CMS simulation)}");
    c->Modified(); c->Update();
    c->SaveAs(Form("../outputDisplay/TriggEff/PairTypeStages_SingleMass__Trigger_%d__%s_bis.pdf", mass, ofilename));

    return;
}

void ProfileVsRunNumber() {

    TFile *ifile = new TFile("/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/JetMET2024_V12/JetMET2024_V12p32.root", "READ");

    // --- Fpix ---
    TH2F *h_Fpix_nosel = (TH2F*)ifile->Get("Nosel_Fpix_vs_RunNumber");
    TH2F *h_Fpix_sel   = (TH2F*)ifile->Get("METanalysis_TestPUppiMETCut_Eta2p4_Fpix_vs_RunNumber");

    // --- Ih ---
    TH2F *h_Ih_nosel = (TH2F*)ifile->Get("Nosel_Ih_vs_RunNumber");
    TH2F *h_Ih_sel   = (TH2F*)ifile->Get("METanalysis_TestPUppiMETCut_Eta2p4_Ih_vs_RunNumber");

    // --- ProfileX (mean of Y per X bin = mean Fpix/Ih per run number) ---
    TProfile *p_Fpix_nosel = h_Fpix_nosel->ProfileX("p_Fpix_nosel");
    TProfile *p_Fpix_sel   = h_Fpix_sel->ProfileX("p_Fpix_sel");
    TProfile *p_Ih_nosel   = h_Ih_nosel->ProfileX("p_Ih_nosel");
    TProfile *p_Ih_sel     = h_Ih_sel->ProfileX("p_Ih_sel");

    p_Fpix_nosel->Rebin(5);
    p_Fpix_sel->Rebin(30);
    p_Ih_nosel->Rebin(5);
    p_Ih_sel->Rebin(30);

    // --- Fill a TH1D with the distribution of the profile bin values ---
    // For each run-number bin, take the profile mean and fill it into a 1D histo.
    // X axis = mean value (Fpix / Ih), Y axis = number of runs.
    auto profileToHist = [](TProfile *p, const char *name, int nbins, double xlo, double xhi) {
        TH1D *h = new TH1D(name, "", nbins, xlo, xhi);
        for (int i = 1; i <= p->GetNbinsX(); ++i) {
            if (p->GetBinEntries(i) > 0)
                h->Fill(p->GetBinContent(i));
        }
        return h;
    };

    TH1D *hp_Fpix_nosel = profileToHist(p_Fpix_nosel, "hp_Fpix_nosel", 50, 0.55, 0.75);
    TH1D *hp_Fpix_sel   = profileToHist(p_Fpix_sel,   "hp_Fpix_sel",   50, 0.55, 0.75);
    TH1D *hp_Ih_nosel   = profileToHist(p_Ih_nosel,   "hp_Ih_nosel",   50, 3.25, 3.5);
    TH1D *hp_Ih_sel     = profileToHist(p_Ih_sel,     "hp_Ih_sel",     50, 3.25, 3.5);

    // --- Styling helper ---
    auto styleProf = [](TProfile *p, int color, int marker, const char *ytitle) {
        p->SetTitle("");
        p->GetXaxis()->SetTitle("Run number");
        p->GetYaxis()->SetTitle(ytitle);
        p->GetYaxis()->SetTitleSize(0.06);
        p->GetXaxis()->SetTitleSize(0.06);
        p->GetXaxis()->SetTitleOffset(0.9);
        p->GetYaxis()->SetTitleOffset(1.0);
        p->GetXaxis()->SetLabelSize(0.04);
        p->GetYaxis()->SetLabelSize(0.05);
        p->SetLineColor(color);
        p->SetMarkerColor(color);
        p->SetMarkerStyle(marker);
        p->SetMarkerSize(0.8);
        p->SetLineWidth(2);
    };

    styleProf(p_Fpix_nosel, kAzure+2, 20, "#LT F_{pixel}#GT");
    styleProf(p_Fpix_sel,   kAzure+2,  21, "#LT F_{pixel}#GT");
    styleProf(p_Ih_nosel,   kAzure+2, 20, "#LT I_{h}#GT [MeV/cm]");
    styleProf(p_Ih_sel,     kAzure+2,  21, "#LT I_{h}#GT [MeV/cm]");

    // --- Styling helper for the distribution histograms ---
    auto styleHist = [](TH1D *h, int color, const char *xtitle) {
        h->SetTitle("");
        h->GetXaxis()->SetTitle(xtitle);
        h->GetYaxis()->SetTitle("Entries");
        h->GetYaxis()->SetTitleSize(0.06);
        h->GetXaxis()->SetTitleSize(0.06);
        h->GetXaxis()->SetTitleOffset(0.9);
        h->GetYaxis()->SetTitleOffset(1);
        h->GetXaxis()->SetLabelSize(0.05);
        h->GetYaxis()->SetLabelSize(0.05);
        h->SetLineColor(color);
        h->SetLineWidth(2);
    };

    styleHist(hp_Fpix_nosel, kAzure+2, "#LT F_{pixel}#GT per run");
    styleHist(hp_Fpix_sel,   kAzure+2,  "#LT F_{pixel}#GT per run");
    styleHist(hp_Ih_nosel,   kAzure+2, "#LT I_{h}#GT per run [MeV/cm]");
    styleHist(hp_Ih_sel,     kAzure+2,  "#LT I_{h}#GT per run [MeV/cm]");

    TLatex *tex = new TLatex(0.68, 0.91, "109 fb^{-1} (13.6 TeV)");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);

    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);

    gStyle->SetOptStat(0);

    // ================= Canvas 1: Fpix with and without selection =================
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    c1->SetLeftMargin(0.16); c1->SetBottomMargin(0.16);

    //p_Fpix_nosel->Draw("E1");
    p_Fpix_sel->Draw("E1");
    p_Fpix_sel->GetYaxis()->SetRangeUser(0.55, 0.75);

    TLegend *leg1 = new TLegend(0.60, 0.72, 0.89, 0.89);
    leg1->SetBorderSize(0);
    leg1->SetTextFont(42);
    leg1->AddEntry(p_Fpix_nosel, "No selection", "lep");
    leg1->AddEntry(p_Fpix_sel,   "HSCP pre-selection", "lep");
    //leg1->Draw();

    tex->Draw();
    latex1->Draw();
    c1->SaveAs("../outputDisplay/ProfileVsRunNumber_Fpix.pdf");
    latex1->SetTitle("#it{Private work (CMS data)}");
    c1->Modified(); c1->Update();
    c1->SaveAs("../outputDisplay/ProfileVsRunNumber_Fpix_bis.pdf");
    latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Work in progress}");

    // ================= Canvas 2: Ih with and without selection =================
    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    c2->SetLeftMargin(0.16); c2->SetBottomMargin(0.16);

    //p_Ih_nosel->Draw("E1");
    p_Ih_sel->Draw("E1");
    p_Ih_sel->GetYaxis()->SetRangeUser(3.25, 3.5);
    

    TLegend *leg2 = new TLegend(0.60, 0.72, 0.89, 0.89);
    leg2->SetBorderSize(0);
    leg2->SetTextFont(42);
    leg2->AddEntry(p_Ih_nosel, "No selection", "lep");
    leg2->AddEntry(p_Ih_sel,   "HSCP pre-selection", "lep");
    //leg2->Draw();

    tex->Draw();
    latex1->Draw();
    c2->SaveAs("../outputDisplay/ProfileVsRunNumber_Ih.pdf");
    latex1->SetTitle("#it{Private work (CMS data)}");
    c2->Modified(); c2->Update();
    c2->SaveAs("../outputDisplay/ProfileVsRunNumber_Ih_bis.pdf");
    latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Work in progress}");

    // ============ Canvas 3: distribution of per-run Fpix profile means ============
    TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
    c3->SetLeftMargin(0.16); c3->SetBottomMargin(0.16);

    //hp_Fpix_nosel->Draw("HIST");
    hp_Fpix_sel->Draw("HIST");
    hp_Fpix_sel->GetYaxis()->SetRangeUser(0, 30);
    

    TLegend *leg3 = new TLegend(0.60, 0.72, 0.89, 0.89);
    leg3->SetBorderSize(0);
    leg3->SetTextFont(42);
    leg3->AddEntry(hp_Fpix_nosel, "No selection", "l");
    leg3->AddEntry(hp_Fpix_sel,   "HSCP pre-selection", "l");
    //leg3->Draw();

    tex->Draw();
    latex1->Draw();
    c3->SaveAs("../outputDisplay/ProfileVsRunNumber_Fpix_hist.pdf");
    latex1->SetTitle("#it{Private work (CMS data)}");
    c3->Modified(); c3->Update();
    c3->SaveAs("../outputDisplay/ProfileVsRunNumber_Fpix_hist_bis.pdf");
    latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Work in progress}");

    // ============ Canvas 4: distribution of per-run Ih profile means ============
    TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);
    c4->SetLeftMargin(0.16); c4->SetBottomMargin(0.16);

    //hp_Ih_nosel->Draw("HIST");
    hp_Ih_sel->Draw("HIST");
    hp_Ih_sel->GetYaxis()->SetRangeUser(0, 30);
    

    TLegend *leg4 = new TLegend(0.60, 0.72, 0.89, 0.89);
    leg4->SetBorderSize(0);
    leg4->SetTextFont(42);
    leg4->AddEntry(hp_Ih_nosel, "No selection", "l");
    leg4->AddEntry(hp_Ih_sel,   "HSCP pre-selection", "l");
    //leg4->Draw();

    tex->Draw();
    latex1->Draw();
    c4->SaveAs("../outputDisplay/ProfileVsRunNumber_Ih_hist.pdf");
    latex1->SetTitle("#it{Private work (CMS data)}");
    c4->Modified(); c4->Update();
    c4->SaveAs("../outputDisplay/ProfileVsRunNumber_Ih_hist_bis.pdf");
    latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Work in progress}");

    return;
}


TCanvas* DrawWithCDF1(TH1* h,
                      TCanvas* c1,
                      std::string CanvasTitle,
                      std::string XaxisTitle,
                      std::string leg_h,
                      float Xmin,
                      float Xmax,
                      int color, bool islogY=false) {

    TCanvas* c_new = new TCanvas(CanvasTitle.c_str(), CanvasTitle.c_str(), 800, 600);

    // Define pads
    TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0);
    pad1->SetLeftMargin(0.16);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();

    TPad* pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.315);
    pad2->SetLeftMargin(0.16); pad2->SetBottomMargin(0.33);
    pad2->Draw();

    // Draw upper plot (the pre-drawn distribution)
    pad1->cd();
    c1->DrawClonePad();
    if (islogY) gPad->SetLogy();

    // Clone histogram
    TH1D* hc = (TH1D*)h->Clone(TString(h->GetName()) + "_" + CanvasTitle.c_str() + "_c");

    // Normalize clone before computing CDF (so CDF goes from 0 to 1)
    if (hc->Integral() > 0) hc->Scale(1.0 / hc->Integral());

    // Build CDF by cumulative sum
    TH1* hCDF = hc->GetCumulative();
    hCDF->SetName(TString("hCDF_") + CanvasTitle.c_str());

    // Style CDF
    hCDF->SetLineColor(color);
    hCDF->SetMarkerColor(color);
    hCDF->SetMarkerStyle(20);

    // Draw CDF in lower pad
    pad2->cd();
    gStyle->SetOptStat(0);
    gPad->SetTickx(0);

    hCDF->SetTitle("");
    hCDF->GetYaxis()->SetTitle("CDF");
    hCDF->GetXaxis()->SetTitle(XaxisTitle.c_str());
    hCDF->GetYaxis()->SetRangeUser(0, 1);
    hCDF->GetYaxis()->SetNdivisions(505);
    hCDF->GetYaxis()->SetTitleFont(43);
    hCDF->GetXaxis()->SetTitleFont(43);
    hCDF->GetYaxis()->SetLabelFont(43);
    hCDF->GetXaxis()->SetLabelFont(43);
    hCDF->GetYaxis()->SetTitleSize(24);
    hCDF->GetXaxis()->SetTitleSize(24);
    hCDF->GetYaxis()->SetLabelSize(20);
    hCDF->GetXaxis()->SetLabelSize(20);
    hCDF->GetYaxis()->SetTitleOffset(1.3);
    hCDF->GetXaxis()->SetTitleOffset(1.0);
    hCDF->LabelsOption("v", "X");
    hCDF->GetXaxis()->SetRangeUser(Xmin, Xmax);

    // draw horizontal dashed line at y=0.5
    TLine* line = new TLine(Xmin, 0.5, Xmax, 0.5);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);

    hCDF->Draw("hist P");
    line->Draw("same");

    c_new->Update();
    cout << "Canvas " << CanvasTitle << " drawn with CDF of: "
         << h->GetName() << endl;

    return c_new;
}

TCanvas* MakeTopCanvas(TH1D *h, const char *cname, const char *legLabel,
                       TLatex *tex, TLatex *latex1, float Xmax, int color) {
    TCanvas *c = new TCanvas(cname, cname, 800, 600);
    c->SetLeftMargin(0.16); c->SetBottomMargin(0.16);

    h->Draw("HIST");
    h->GetYaxis()->SetRangeUser(0.1, 1.3 * h->GetMaximum());
    h->GetXaxis()->SetRangeUser(0, Xmax);
    h->SetLineColor(color);

    tex->Draw();
    latex1->Draw();
    c->Modified(); c->Update();
    return c;
}

void CompareIhWithCDF() {

    TFile *ifileData   = new TFile("../output/JetMET2024_V12/JetMET2024_V12p31.root", "READ");
    TFile *ifileSignal = new TFile("../output/Gluino_V19/Gluino_Run3_MET_madgraph_2000_V19p7.root", "READ");

    TH1D *IhSignal   = (TH1D*)ifileSignal->Get("METanalysis_TestPUppiMETCut_Eta2p4_Ih");
    TH1D *IhDataRaw  = (TH1D*)ifileData->Get("Nm1_Ih_StripOnly");

    // --- Rebin data (200 bins, 0-10) onto signal binning (600 bins, 0-30) ---
    TH1D *IhData = (TH1D*)IhSignal->Clone("IhData");
    IhData->Reset();
    for (int i = 1; i <= IhDataRaw->GetNbinsX(); ++i) {
        double x = IhDataRaw->GetBinCenter(i);
        int j = IhData->FindBin(x);
        IhData->AddBinContent(j, IhDataRaw->GetBinContent(i));
    }
    for (int j = 1; j <= IhData->GetNbinsX(); ++j)
        IhData->SetBinError(j, TMath::Sqrt(IhData->GetBinContent(j)));

    // --- Styling ---
    auto styleHist = [](TH1D *h, int color, int marker) {
        h->SetTitle("");
        h->GetXaxis()->SetTitle("I_{h} [MeV/cm]");
        h->GetYaxis()->SetTitle("Normalised entries");
        h->GetYaxis()->SetTitleSize(0.06);
        h->GetXaxis()->SetTitleSize(0.06);
        h->GetXaxis()->SetTitleOffset(0.9);
        h->GetYaxis()->SetTitleOffset(1.0);
        h->GetXaxis()->SetLabelSize(0.05);
        h->GetYaxis()->SetLabelSize(0.05);
        h->SetLineColor(color);
        h->SetMarkerColor(color);
        h->SetMarkerStyle(marker);
        h->SetMarkerSize(0.8);
        h->SetLineWidth(2);
    };

    styleHist(IhData,   kAzure+2, 20);
    styleHist(IhSignal, kRed+1,   21);

    TLatex *tex = new TLatex(0.68, 0.91, "109 fb^{-1} (13.6 TeV)");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);

    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{Private work (CMS data)}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);

    gStyle->SetOptStat(0);

    // ================= Canvas 1: Data + its CDF =================
    TCanvas *cTopData = MakeTopCanvas(IhData, "cTopData", "Data", tex, latex1, 7., 602);
    TCanvas *cData = DrawWithCDF1(IhData, cTopData,
                                  "Ih_Data_CDF",
                                  "I_{h} [MeV/cm]",
                                  "Data",
                                  0.0, 7.,
                                  602, true);
    cData->SaveAs("../outputDisplay/Ih_Data_CDF.pdf");

    latex1->SetTitle("#it{Private work (CMS simulation)}");
    tex->SetTitle("");

    // ================= Canvas 2: Signal + its CDF =================
    TCanvas *cTopSignal = MakeTopCanvas(IhSignal, "cTopSignal", "Gluino 2000 GeV", tex, latex1, 30, kRed+1);
    TCanvas *cSignal = DrawWithCDF1(IhSignal, cTopSignal,
                                    "Ih_Signal_CDF",
                                    "I_{h} [MeV/cm]",
                                    "Gluino 2000 GeV",
                                    0.0, 30.0,
                                    kRed+1, true);
    cSignal->SaveAs("../outputDisplay/Ih_Signal_CDF.pdf");

    return;
}


void PUppiVSPF(std::string ifileName, bool ifData) {

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);

    TFile *ifile = new TFile(ifileName.c_str(), "READ");
    if (!ifile || ifile->IsZombie()) {
        std::cerr << "Error: cannot open " << ifileName << std::endl;
        return;
    }

    TH2F *h2[2];
    h2[0] = (TH2F*)ifile->Get("Nosel_PUppiMET_VS_PseudoMET");          // nosel
    h2[1] = (TH2F*)ifile->Get("HSCPPartialsel_PUppiMET_VS_PseudoMET"); // sel
    if (!h2[0] || !h2[1]) {
        std::cerr << "Error: missing 2D histos in " << ifileName << std::endl;
        return;
    }

    // ================= TLatex (habillage de base) =================
    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    if (ifData) latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);

    TLatex *latex2 = new TLatex(0.68, 0.91, "#scale[1.3]{#bf{m_{#tilde{g}}=2000 GeV}}");
    if (ifData) latex2->SetTitle("109 fb^{-1} (13.6 TeV)");
    latex2->SetNDC();
    latex2->SetTextFont(42);
    latex2->SetTextSize(0.04);

    std::cout << "nosel entries = " << h2[0]->GetEntries()
          << "  integral = " << h2[0]->Integral() << std::endl;
    std::cout << "sel   entries = " << h2[1]->GetEntries()
          << "  integral = " << h2[1]->Integral() << std::endl;

    // ================= Canvas 2D COLZ (sel + nosel) =================
    for (int s = 0; s < 2; ++s) {
        const char *tag = (s == 0) ? "nosel" : "sel";

        h2[s]->SetTitle("");
        h2[s]->GetYaxis()->SetTitle("Pseudo MET [GeV]");
        h2[s]->GetXaxis()->SetTitle("PUppi MET [GeV]");
        h2[s]->GetZaxis()->SetTitle("Entries");
        h2[s]->GetYaxis()->SetTitleSize(0.06);
        h2[s]->GetXaxis()->SetTitleSize(0.06);
        h2[s]->GetXaxis()->SetTitleOffset(0.9);
        h2[s]->GetYaxis()->SetTitleOffset(1.1);
        h2[s]->GetXaxis()->SetLabelSize(0.05);
        h2[s]->GetYaxis()->SetLabelSize(0.05);
        h2[s]->GetZaxis()->SetTitleSize(0.06);
        h2[s]->GetZaxis()->SetTitleOffset(0.9);
        h2[s]->GetZaxis()->SetLabelSize(0.05);
        if (!ifData) {
            h2[s]->SetMinimum(0.001);
            h2[s]->SetMaximum(5);
        }

        TCanvas *c6 = new TCanvas(Form("c6_%s", tag), "c6", 800, 600);
        c6->cd();
        c6->SetLeftMargin(0.16); c6->SetBottomMargin(0.16); c6->SetRightMargin(0.16);

        h2[s]->Draw("COLZ");
        c6->SetLogz();
        latex1->Draw();
        latex2->Draw();
        c6->SaveAs(Form("../outputDisplay/PUppiVSPF_%s_%s.pdf", (ifData ? "Data" : "MC"), tag));

        TString saved = latex1->GetTitle();
        latex1->SetTitle((ifData) ? "#it{Private work (CMS data)}" : "#it{Private work (CMS simulation)}");
        c6->Modified(); c6->Update();
        c6->SaveAs(Form("../outputDisplay/PUppiVSPF_%s_%s_bis.pdf", (ifData ? "Data" : "MC"), tag));
        latex1->SetTitle(saved);
    }

    // ================= Projections 1D =================
    // [0] = nosel (trait plein, marker plein), [1] = sel (pointille, marker vide)
    TH1D *hPUppi[2], *hPseudo[2];
    hPUppi[0]  = h2[0]->ProjectionX("PUppi_nosel");
    hPseudo[0] = h2[0]->ProjectionY("Pseudo_nosel");
    hPUppi[1]  = h2[1]->ProjectionX("PUppi_sel");
    hPseudo[1] = h2[1]->ProjectionY("Pseudo_sel");
    for (int s = 0; s < 2; ++s) { hPUppi[s]->SetDirectory(0); hPseudo[s]->SetDirectory(0); }

    const int colPUppi  = kOrange+8;
    const int colPseudo = kGreen+3;

    for (int s = 0; s < 2; ++s) {
        int lstyle = (s == 0) ? 1 : 2;
        hPUppi[s]->SetLineColor(colPUppi);   hPUppi[s]->SetMarkerColor(colPUppi);
        hPUppi[s]->SetLineWidth(2);          hPUppi[s]->SetLineStyle(lstyle);
        hPseudo[s]->SetLineColor(colPseudo); hPseudo[s]->SetMarkerColor(colPseudo);
        hPseudo[s]->SetLineWidth(2);         hPseudo[s]->SetLineStyle(lstyle);
    }

    // habillage de base sur l'histo porteur du pad haut
    hPUppi[0]->SetTitle("");
    hPUppi[0]->GetYaxis()->SetTitle("Entries");
    hPUppi[0]->GetXaxis()->SetTitle("MET [GeV]");
    hPUppi[0]->GetYaxis()->SetTitleSize(0.06);
    hPUppi[0]->GetYaxis()->SetTitleOffset(1.1);
    hPUppi[0]->GetYaxis()->SetLabelSize(0.05);
    hPUppi[0]->GetXaxis()->SetLabelSize(0.05);   // axe X porte par le pad du bas
    hPUppi[0]->GetXaxis()->SetTitleSize(0.06);
    hPUppi[0]->GetXaxis()->SetTitleOffset(1.0);

    double ymax = 0.;
    for (int s = 0; s < 2; ++s) {
        ymax = std::max(ymax, hPUppi[s]->GetMaximum());
        ymax = std::max(ymax, hPseudo[s]->GetMaximum());
    }
    hPUppi[0]->SetMaximum(2 * ymax);
    hPUppi[0]->SetMinimum(0.1);

    // ================= CDF (4 courbes) =================
    // markers : nosel plein (20/21), sel vide (24/25)
    auto makeCDF = [](TH1D *h, int color, int mstyle, int lstyle) {
        TH1D *hc = (TH1D*)h->Clone(Form("%s_cdfclone", h->GetName()));
        hc->SetDirectory(0);
        if (hc->Integral() > 0) hc->Scale(1.0 / hc->Integral());
        TH1 *cdf = hc->GetCumulative();
        cdf->SetName(Form("%s_cdf", h->GetName()));
        cdf->SetDirectory(0);
        cdf->SetLineColor(color);    cdf->SetMarkerColor(color);
        cdf->SetMarkerStyle(mstyle); cdf->SetMarkerSize(0.8);
        cdf->SetLineWidth(2);        cdf->SetLineStyle(lstyle);
        return cdf;
    };

    TH1 *cdfPUppi[2], *cdfPseudo[2];
    cdfPUppi[0]  = makeCDF(hPUppi[0],  colPUppi,  20, 1);  // nosel plein
    cdfPseudo[0] = makeCDF(hPseudo[0], colPseudo, 21, 1);
    cdfPUppi[1]  = makeCDF(hPUppi[1],  colPUppi,  24, 2);  // sel vide
    cdfPseudo[1] = makeCDF(hPseudo[1], colPseudo, 25, 2);

    double Xmin = hPUppi[0]->GetXaxis()->GetXmin();
    double Xmax = hPUppi[0]->GetXaxis()->GetXmax();

    // ================= Canvas projections + CDF (2 pads) =================
    TCanvas *c7 = new TCanvas("c7", "c7", 800, 600);

    TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0);
    pad1->SetLeftMargin(0.16);
    pad1->SetRightMargin(0.16);
    pad1->SetBottomMargin(0.15);
    pad1->Draw();

    TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.315);
    pad2->SetLeftMargin(0.16);
    pad2->SetRightMargin(0.16);
    pad2->SetBottomMargin(0.33);
    pad2->SetTopMargin(0.03);
    pad2->Draw();

    // --- pad du haut : distributions ---
    pad1->cd();
    pad1->SetLogy();

    hPUppi[0]->Draw("HIST");
    hPUppi[1]->Draw("HIST SAME");
    hPseudo[0]->Draw("HIST SAME");
    hPseudo[1]->Draw("HIST SAME");

    TLegend *leg = new TLegend(0.55, 0.60, 0.85, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hPUppi[0],  "PUppi MET (no sel)",  "l");
    leg->AddEntry(hPUppi[1],  "PUppi MET (sel)",     "l");
    leg->AddEntry(hPseudo[0], "Pseudo MET (no sel)", "l");
    leg->AddEntry(hPseudo[1], "Pseudo MET (sel)",    "l");
    leg->Draw();

    latex1->Draw();
    latex2->Draw();

    // --- pad du bas : CDF ---
    pad2->cd();
    gPad->SetTickx(0);

    cdfPUppi[0]->SetTitle("");
    cdfPUppi[0]->GetYaxis()->SetTitle("CDF");
    cdfPUppi[0]->GetXaxis()->SetTitle("MET [GeV]");
    cdfPUppi[0]->GetYaxis()->SetRangeUser(0, 1);
    cdfPUppi[0]->GetYaxis()->SetNdivisions(505);
    cdfPUppi[0]->GetYaxis()->SetTitleFont(43);
    cdfPUppi[0]->GetXaxis()->SetTitleFont(43);
    cdfPUppi[0]->GetYaxis()->SetLabelFont(43);
    cdfPUppi[0]->GetXaxis()->SetLabelFont(43);
    cdfPUppi[0]->GetYaxis()->SetTitleSize(24);
    cdfPUppi[0]->GetXaxis()->SetTitleSize(24);
    cdfPUppi[0]->GetYaxis()->SetLabelSize(20);
    cdfPUppi[0]->GetXaxis()->SetLabelSize(20);
    cdfPUppi[0]->GetYaxis()->SetTitleOffset(1.3);
    cdfPUppi[0]->GetXaxis()->SetTitleOffset(1.0);
    cdfPUppi[0]->GetXaxis()->SetRangeUser(Xmin, Xmax);

    cdfPUppi[0]->Draw("hist P");
    cdfPUppi[1]->Draw("hist P SAME");
    cdfPseudo[0]->Draw("hist P SAME");
    cdfPseudo[1]->Draw("hist P SAME");

    TLine *line = new TLine(Xmin, 0.5, Xmax, 0.5);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->Draw("same");

    // ================= Sauvegardes =================
    c7->cd();
    c7->Modified(); c7->Update();
    c7->SaveAs(Form("../outputDisplay/PUppiPseudoProj_%s.pdf", (ifData ? "Data" : "MC")));

    latex1->SetTitle((ifData) ? "#it{Private work (CMS data)}" : "#it{Private work (CMS simulation)}");
    c7->Modified(); c7->Update();
    c7->SaveAs(Form("../outputDisplay/PUppiPseudoProj_%s_bis.pdf", (ifData ? "Data" : "MC")));

    return;
}


TCanvas* DrawWithCDF(TH1F* data, TH1F* mc, TH1F* sig, TCanvas* cMain,
                     TString cname, TString xlabel, double xmin, double xmax) {

    TCanvas* c = new TCanvas(cname, cname, 800, 600);

    // --- pad du haut : on clone le plot déjà construit ---
    TPad* pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetLeftMargin(0.16);
    pad1->Draw();
    pad1->cd();
    cMain->SetBottomMargin(0.13);
    cMain->Modified();
    cMain->DrawClonePad();

    // FORCE le range X sur tous les TH1 clonés (frame inclus)
    TIter next(pad1->GetListOfPrimitives());
    TObject* obj;
    while ((obj = next())) {
        if (obj->InheritsFrom(TH1::Class()))
            ((TH1*)obj)->GetXaxis()->SetRangeUser(xmin, xmax);
    }
    pad1->Modified();
    pad1->Update();


    // --- pad du bas : CDF ---
    c->cd();
    TPad* pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);
    pad2->SetLeftMargin(0.16);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.33);
    pad2->SetTickx(1);
    pad2->SetTicky(1);
    pad2->Draw();
    pad2->cd();

    // CDF = cumulative normalisée à 1
    TH1F* cdfData = (TH1F*)data->GetCumulative();
    TH1F* cdfMC   = (TH1F*)mc->GetCumulative();
    TH1F* cdfSig  = (TH1F*)sig->GetCumulative();

    double iData = data->Integral(1, data->GetNbinsX());
    double iMC   = mc->Integral(1, mc->GetNbinsX());
    double iSig  = sig->Integral(1, sig->GetNbinsX());
    if (iData > 0) cdfData->Scale(1.0/iData);
    if (iMC   > 0) cdfMC->Scale(1.0/iMC);
    if (iSig  > 0) cdfSig->Scale(1.0/iSig);

    // style
    cdfMC->SetLineColor(kBlue-7);
    cdfMC->SetLineWidth(2);
    cdfMC->SetMarkerStyle(21);
    cdfMC->SetMarkerColor(kBlue-7);

    cdfSig->SetLineColor(kViolet-1);
    cdfSig->SetLineWidth(2);
    cdfSig->SetMarkerStyle(22);
    cdfSig->SetMarkerColor(kViolet-1);

    cdfData->SetLineColor(kBlack);
    cdfData->SetLineWidth(2);
    cdfData->SetMarkerStyle(20);
    cdfData->SetMarkerColor(kBlack);

    cdfMC->SetTitle("");
    cdfMC->SetMinimum(0.0);
    cdfMC->SetMaximum(1.1);
    cdfMC->GetXaxis()->SetRangeUser(xmin, xmax);
    cdfData->GetXaxis()->SetRangeUser(xmin, xmax);
    cdfSig->GetXaxis()->SetRangeUser(xmin, xmax);
    cdfMC->GetXaxis()->SetTitle(xlabel);
    cdfMC->GetYaxis()->SetTitle("CDF");
    cdfMC->GetYaxis()->SetNdivisions(505);
    cdfMC->GetXaxis()->SetTitleSize(0.15);
    cdfMC->GetXaxis()->SetTitleOffset(0.9);
    cdfMC->GetXaxis()->SetLabelSize(0.12);
    cdfMC->GetYaxis()->SetTitleSize(0.14);
    cdfMC->GetYaxis()->SetTitleOffset(0.3);
    cdfMC->GetYaxis()->SetLabelSize(0.12);

    cdfMC->Draw("E0");
    cdfSig->Draw("E0 same");
    cdfData->Draw("E0 same");

    TLine* line = new TLine(xmin, 0.5, xmax, 0.5);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->Draw("same");

    pad2->RedrawAxis();
    c->cd();
    c->Update();
    return c;
}


void Nm1Eff (bool isRescaled,
             std::string TTbar,
             std::string TTbarSemiLep,
             std::string Wjets,
             std::string QCD,
             std::string JetMETdata,
             std::string Gluino) {

    gErrorIgnoreLevel = kWarning;

    TFile *ifile_TTbar = new TFile(TTbar.c_str(), "READ");
    TFile *ifile_TTbarSemiLep = new TFile(TTbarSemiLep.c_str(), "READ");
    TFile *ifile_Wjets = new TFile(Wjets.c_str(), "READ");
    TFile *ifile_QCD = new TFile(QCD.c_str(), "READ");
    TFile *ifile_JetMETdata = new TFile(JetMETdata.c_str(), "READ");
    TFile *ifile_Gluino = new TFile(Gluino.c_str(), "READ");

    std::vector<TH1F*> hNM1_TTbar, hNM1_TTbarSemiLep, hNM1_Wjets, hNM1_QCD, hNM1_MC, hNM1_JetMETdata, hNM1_Gluino;

    std::vector<string> cutNames = {"trigger", "METfilters", "PUppiMET", "Ptpseudo", "eta", "NOPH", "NOM", "FOVH", "HighPurity", "Chi2", "dZ", "dXY", "PFMiniIso", "TrkIso", "EoverP", "PtErr_over_PtPt", "PtErr_over_Pt", "Fpix", "Ih_StripOnly", "Ih_StripOnly_rescaled"};

    std::vector<string> Xlabel = {"or MET trigger", "MET filters", "PUppi MET", "p_{T} [GeV]", "#eta", "Nb pixel hits", "Nb dE/dx", "frac. valid hits", "HighPurity", "#chi^{2}/NDOF", "d_{z} [cm]", "d_{xy} [cm]", "I_{PF}^{rel}", "I_{dr03}^{trk}", "E/p", "#sigma_{p_{T}}/p_{T}^{2}", "#sigma_{p_{T}}/p_{T}", "F_{pixel}", "I_{h} [MeV/cm]", "I_{h} rescaled [MeV/cm]"};

    // Step 1: retrieve the Nm1 hists
    for (size_t i = 0; i < cutNames.size(); i++) {

        if (cutNames[i] == "Ih_StripOnly_rescaled") {
            hNM1_JetMETdata.push_back((TH1F*)ifile_JetMETdata->Get(Form("Nm1_event_%s", cutNames[i-1].c_str())));
            hNM1_TTbar.push_back((TH1F*)ifile_TTbar->Get(Form("Nm1_event_%s", cutNames[i].c_str())));
            hNM1_TTbarSemiLep.push_back((TH1F*)ifile_TTbarSemiLep->Get(Form("Nm1_event_%s", cutNames[i].c_str())));
            hNM1_Wjets.push_back((TH1F*)ifile_Wjets->Get(Form("Nm1_event_%s", cutNames[i].c_str())));
            hNM1_QCD.push_back((TH1F*)ifile_QCD->Get(Form("Nm1_event_%s", cutNames[i].c_str())));
            hNM1_Gluino.push_back((TH1F*)ifile_Gluino->Get(Form("Nm1_event_%s", cutNames[i].c_str())));
            hNM1_MC.push_back((TH1F*)ifile_Wjets->Get(Form("Nm1_event_%s", cutNames[i].c_str())));
        }

        hNM1_JetMETdata.push_back((TH1F*)ifile_JetMETdata->Get(Form("Nm1_event_%s", cutNames[i].c_str())));
        hNM1_TTbar.push_back((TH1F*)ifile_TTbar->Get(Form("Nm1_event_%s", cutNames[i].c_str())));
        hNM1_TTbarSemiLep.push_back((TH1F*)ifile_TTbarSemiLep->Get(Form("Nm1_event_%s", cutNames[i].c_str())));
        hNM1_Wjets.push_back((TH1F*)ifile_Wjets->Get(Form("Nm1_event_%s", cutNames[i].c_str())));
        hNM1_QCD.push_back((TH1F*)ifile_QCD->Get(Form("Nm1_event_%s", cutNames[i].c_str())));
        hNM1_Gluino.push_back((TH1F*)ifile_Gluino->Get(Form("Nm1_event_%s", cutNames[i].c_str())));
        hNM1_MC.push_back((TH1F*)ifile_Wjets->Get(Form("Nm1_event_%s", cutNames[i].c_str())));
    }

    // Step 2: style + build MC sum
    for (size_t i = 0; i < cutNames.size(); i++) {
        cout << "Cut: " << cutNames[i] << endl;
        hNM1_MC[i]->Add(hNM1_TTbar[i]);
        hNM1_MC[i]->Add(hNM1_QCD[i]);

        if (i==2 || i==3 || i==13 || i==14 || i==15 || i==16) {
            hNM1_TTbar[i]->Rebin(4);
            hNM1_TTbarSemiLep[i]->Rebin(4);
            hNM1_Wjets[i]->Rebin(4);
            hNM1_QCD[i]->Rebin(4);
            hNM1_Gluino[i]->Rebin(4);
            hNM1_JetMETdata[i]->Rebin(4);
        }
        if (i==12) {
            hNM1_TTbar[i]->Rebin(8);
            hNM1_TTbarSemiLep[i]->Rebin(8);
            hNM1_Wjets[i]->Rebin(8);
            hNM1_QCD[i]->Rebin(8);
            hNM1_Gluino[i]->Rebin(8);
            hNM1_JetMETdata[i]->Rebin(8);
        }

        hNM1_TTbar[i]->SetLineColor(kRed);
        hNM1_TTbar[i]->SetFillColorAlpha(kRed, 0.5);
        hNM1_TTbarSemiLep[i]->SetLineColor(kRed+2);
        hNM1_TTbarSemiLep[i]->SetFillColorAlpha(kRed+2, 0.5);
        hNM1_Wjets[i]->SetLineColor(kBlue-7);
        hNM1_Wjets[i]->SetFillColorAlpha(kBlue-7, 0.5);
        hNM1_QCD[i]->SetLineColor(kGreen);
        hNM1_QCD[i]->SetFillColorAlpha(kGreen, 0.5);
        hNM1_Gluino[i]->SetLineColor(kViolet-1);
        hNM1_Gluino[i]->SetMarkerStyle(22);
        hNM1_Gluino[i]->SetMarkerColor(kViolet-1);
        hNM1_JetMETdata[i]->SetLineColor(kBlack);
        hNM1_JetMETdata[i]->SetMarkerColor(kBlack);
        hNM1_JetMETdata[i]->SetMarkerStyle(20);
    }

    // Step 2 bis: if rescaled, scale MC to data
    if (isRescaled) {
        for (size_t i = 0; i < cutNames.size(); i++) {
            float intTot = hNM1_TTbar[i]->Integral() + hNM1_TTbarSemiLep[i]->Integral() + hNM1_Wjets[i]->Integral() + hNM1_QCD[i]->Integral();
            hNM1_TTbar[i]->Scale(hNM1_JetMETdata[i]->Integral()/intTot);
            hNM1_TTbarSemiLep[i]->Scale(hNM1_JetMETdata[i]->Integral()/intTot);
            hNM1_Wjets[i]->Scale(hNM1_JetMETdata[i]->Integral()/intTot);
            hNM1_QCD[i]->Scale(hNM1_JetMETdata[i]->Integral()/intTot);
            hNM1_Gluino[i]->Scale(hNM1_JetMETdata[i]->Integral()/hNM1_Gluino[i]->Integral());
        }
    }

    struct CutInfo {
        double threshold;
        bool keepRight;
        bool isSymmetric;
    };

    std::vector<CutInfo> cutInfos = {
        {0.5,      true,  false},  // trigger
        {0.5,      true,  false},  // METfilters
        {150,    true,  false},  // PUppiMET > 150
        {50,     true,  false},  // Pt_pseudo     > 50
        {2.4,    false, true },  // |eta|         < 2.4
        {2,      true,  false},  // NOPH          >= 2
        {10,     true,  false},  // NOM           >= 10
        {0.8,    true,  false},  // FOVH          > 0.8
        {0.5,      true,  false},  // HighPurity
        {5.0,    false, false},  // Chi2          < 5
        {0.1,    false, true },  // |dZ|          < 0.1
        {0.02,   false, true },  // |dXY|         < 0.02
        {0.02,   false, false},  // PFMiniIso     < 0.02
        {15,     false, false},  // TrkIso        < 15
        {0.3,    false, false},  // EoverP        < 0.3
        {0.0008, false, false},  // PtErr/PtPt    < 0.0008
        {1,      false, false},  // PtErr/Pt      < 1
        {0.3,    true,  false},  // Fpix          > 0.3
        {2.9784, true,  false},  // Ih_StripOnly
        {2.9784, true,  false},  // Ih_StripOnly rescaled
    };

    TLatex *tex = new TLatex(0.75, 0.91, "109 fb^{-1} (13.6 TeV)");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);

    TLatex *latex1 = new TLatex(0.16, 0.91, "#it{Private work (CMS simulation/data)}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);

    auto addOverflow = [](TH1F* h) {
        int n = h->GetNbinsX();
        h->SetBinContent(n, h->GetBinContent(n) + h->GetBinContent(n+1));
        h->SetBinError(n, std::sqrt(std::pow(h->GetBinError(n),2) + std::pow(h->GetBinError(n+1),2)));
        h->SetBinContent(n+1, 0);   // évite le double comptage
        h->SetBinError(n+1, 0);
    };

    // Step 3: Canvas for each hists
    std::vector<TCanvas*> canvases;
    for (size_t i = 0; i < cutNames.size(); i++) {

        TCanvas *c = new TCanvas(Form("c_Nm1_%s", cutNames[i].c_str()), Form("c_Nm1_%s", cutNames[i].c_str()), 800, 600);
        c->SetLeftMargin(0.16); c->SetBottomMargin(0.16);

        addOverflow(hNM1_QCD[i]);
        addOverflow(hNM1_TTbar[i]);
        addOverflow(hNM1_TTbarSemiLep[i]);
        addOverflow(hNM1_Wjets[i]);
        addOverflow(hNM1_JetMETdata[i]);
        addOverflow(hNM1_Gluino[i]);

        THStack *hs = new THStack(Form("hs_%s", cutNames[i].c_str()), "");
        hs->Add(hNM1_QCD[i]);
        hs->Add(hNM1_TTbar[i]);
        hs->Add(hNM1_TTbarSemiLep[i]);
        hs->Add(hNM1_Wjets[i]);

        double ymin = 0.01;
        double ymax = 6e6;
        double xmin = hNM1_TTbar[i]->GetBinLowEdge(1);
        double xmax = hNM1_TTbar[i]->GetBinLowEdge(hNM1_TTbar[i]->GetNbinsX()+1);

        if (i==6) xmin = 0.75;
        
        c->cd();
        TH1F *frame = (TH1F*)hNM1_TTbar[i]->Clone(Form("frame_%s", cutNames[i].c_str()));
        gStyle->SetOptStat(0);
        frame->Reset();
        frame->SetMinimum(ymin);
        frame->SetMaximum(ymax);
        frame->GetXaxis()->SetTitle(Xlabel[i].c_str());
        frame->GetYaxis()->SetTitle("Number of events");
        frame->GetYaxis()->SetTitleSize(0.06);
        frame->GetXaxis()->SetTitleSize(0.06);
        frame->GetXaxis()->SetTitleOffset(0.9);
        frame->GetYaxis()->SetTitleOffset(0.8);
        frame->GetXaxis()->SetLabelSize(0.05);
        frame->GetYaxis()->SetLabelSize(0.05);
        frame->GetXaxis()->SetRangeUser(xmin, xmax);

        frame->Draw("HIST");
        hs->Draw("HIST same");
        hNM1_JetMETdata[i]->Draw("E1 same");
        hNM1_Gluino[i]->Draw("E1 same");


        // ligne verticale à la valeur de la coupure
        double xcut = cutInfos[i].threshold;
        TLine *line = new TLine(xcut, ymin, xcut, ymax);
        line->SetLineColor(kBlack);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");

        // pour les coupures symétriques |var| < seuil, tracer aussi la ligne à -seuil
        if (cutInfos[i].isSymmetric) {
            TLine *lineNeg = new TLine(-xcut, ymin, -xcut, ymax);
            lineNeg->SetLineColor(kBlack);
            lineNeg->SetLineStyle(2);
            lineNeg->SetLineWidth(2);
            lineNeg->Draw("same");
        }

        tex->Draw();
        latex1->Draw();


        TLegend *legend = new TLegend(0.7, 0.6, 0.89, 0.89);
        if (i==0 || i==1 || i==6 || i==4 || i==8 || i==17) legend = new TLegend(0.67, 0.16, 0.89, 0.45);
        legend->AddEntry(hNM1_TTbar[i],     "t#bar{t}#rightarrow2l2#nu",     "f");
        legend->AddEntry(hNM1_TTbarSemiLep[i],     "t#bar{t}#rightarrowl#nu",     "f");
        legend->AddEntry(hNM1_Wjets[i],     "W(#rightarrow#mu#nu)+jets",    "f");
        legend->AddEntry(hNM1_QCD[i],       "QCD (#mu enriched)",       "f");
        legend->AddEntry(hNM1_JetMETdata[i],"MET data","lep");
        legend->AddEntry(hNM1_Gluino[i],     "#tilde{g} (m=2000 GeV)", "lep");
        legend->SetBorderSize(0);

        TH1F *htemp = (TH1F*)hNM1_Wjets[i]->Clone("htemp");
        htemp->Add(hNM1_TTbar[i]);
        htemp->Add(hNM1_TTbarSemiLep[i]);
        htemp->Add(hNM1_QCD[i]);

        htemp->SetMinimum(ymin);
        htemp->SetMaximum(ymax);
        hNM1_JetMETdata[i]->SetMinimum(ymin);
        hNM1_JetMETdata[i]->SetMaximum(ymax);

        c->cd();
        gPad->RedrawAxis();
        c->SetLogy();
        c->SetTickx(1);
        c->SetTicky(1);
        hs->SetMinimum(ymin);
        hs->SetMaximum(ymax);
        c->Modified();
        legend->Draw();
        c->Update();

        TCanvas *cRatio = DrawWithCDF(hNM1_JetMETdata[i], htemp, hNM1_Gluino[i], c,
                                      Form("cCDF_%s", cutNames[i].c_str()), Xlabel[i], xmin, xmax);

        canvases.push_back(cRatio);
        delete htemp;
    }

    // Step 4: compute the efficiency for each cut
    auto computeEff = [](TH1F* h, double threshold, bool keepRight, bool isSymmetric = false) -> double {
        if (!h) return -1.0;
        int totalBins = h->GetNbinsX();
        double total  = h->Integral(1, totalBins);
        if (total == 0) return 0.0;
        if (isSymmetric) {
            int binLow  = h->FindBin(-threshold);
            int binHigh = h->FindBin(threshold);
            return h->Integral(binLow, binHigh) / total;
        }
        int cutBin = h->FindBin(threshold);
        double signal = keepRight ? h->Integral(cutBin, totalBins)
                                  : h->Integral(1, cutBin - 1);
        return signal / total;
    };

    // Print efficiencies
    std::cout << std::left
            << std::setw(25) << "Cut"
            << std::setw(12) << "TTbar"
            << std::setw(12) << "TTbarSemiLep"
            << std::setw(12) << "W+jets"
            << std::setw(12) << "QCD"
            << std::setw(12) << "ALL MC"
            << std::setw(12) << "Data"
            << std::endl;
    std::cout << std::string(60, '-') << std::endl;

    for (size_t i = 0; i < cutNames.size(); i++) {
        double eff_TTbar = computeEff(hNM1_TTbar[i],      cutInfos[i].threshold, cutInfos[i].keepRight);
        double eff_TTbarSemiLep = computeEff(hNM1_TTbarSemiLep[i],      cutInfos[i].threshold, cutInfos[i].keepRight);
        double eff_Wjets = computeEff(hNM1_Wjets[i],      cutInfos[i].threshold, cutInfos[i].keepRight);
        double eff_QCD   = computeEff(hNM1_QCD[i],        cutInfos[i].threshold, cutInfos[i].keepRight);
        double eff_MC    = computeEff(hNM1_MC[i],         cutInfos[i].threshold, cutInfos[i].keepRight);
        double eff_Data  = computeEff(hNM1_JetMETdata[i], cutInfos[i].threshold, cutInfos[i].keepRight);

        std::cout << std::left  << std::setw(25) << cutNames[i]
                << std::fixed << std::setprecision(4)
                << std::setw(12) << eff_TTbar
                << std::setw(12) << eff_TTbarSemiLep
                << std::setw(12) << eff_Wjets
                << std::setw(12) << eff_QCD
                << std::setw(12) << eff_MC
                << std::setw(12) << eff_Data
                << std::endl;
    }

    // Step 5: save
    TString pdfName = Form("../outputDisplay/Nm1plots/Nm1Eff%s_ALL.pdf", isRescaled ? "_rescaled" : "");
    for (size_t i = 0; i < canvases.size(); i++) {
        if      (i == 0)                    canvases[i]->Print(pdfName + "(");
        else if (i == canvases.size() - 1)  canvases[i]->Print(pdfName + ")");
        else                                canvases[i]->Print(pdfName);
    }

    return;
}

void NoselEff (bool isRescaled,
               std::string TTbar,
               std::string TTbarSemiLep,
               std::string Wjets,
               std::string QCD,
               std::string JetMETdata,
               std::string Gluino) {

    gErrorIgnoreLevel = kWarning;

    TFile *ifile_TTbar = new TFile(TTbar.c_str(), "READ");
    TFile *ifile_TTbarSemiLep = new TFile(TTbarSemiLep.c_str(), "READ");
    TFile *ifile_Wjets = new TFile(Wjets.c_str(), "READ");
    TFile *ifile_QCD = new TFile(QCD.c_str(), "READ");
    TFile *ifile_JetMETdata = new TFile(JetMETdata.c_str(), "READ");
    TFile *ifile_Gluino = new TFile(Gluino.c_str(), "READ");

    std::vector<TH1F*> hNosel_TTbar, hNosel_TTbarSemiLep, hNosel_Wjets, hNosel_QCD, hNosel_MC, hNosel_JetMETdata, hNosel_Gluino;

    std::vector<string> cutNames = {"PUppiMET", "Ptpseudo", "eta", "PFMiniIso", "TrkIso", "EoverP",
    "PtErr_over_PtPt", "PtErr_over_Pt", "Fpix", "Ih"};

    std::vector<string> Xlabel = {"PUppi MET", "p_{T} [GeV]", "#eta", "I_{PF}^{rel}", "I_{dr03}^{trk}", "E/p",
    "#sigma_{p_{T}}/p_{T}^{2}", "#sigma_{p_{T}}/p_{T}", "F_{pixel}", "I_{h} [MeV/cm]"};

    // Step 1: retrieve the Nosel hists
    for (size_t i = 0; i < cutNames.size(); i++) {
        hNosel_JetMETdata.push_back((TH1F*)ifile_JetMETdata->Get(Form("Nosel_%s", cutNames[i].c_str())));
        hNosel_TTbar.push_back((TH1F*)ifile_TTbar->Get(Form("Nosel_%s", cutNames[i].c_str())));
        hNosel_TTbarSemiLep.push_back((TH1F*)ifile_TTbarSemiLep->Get(Form("Nosel_%s", cutNames[i].c_str())));
        hNosel_Wjets.push_back((TH1F*)ifile_Wjets->Get(Form("Nosel_%s", cutNames[i].c_str())));
        hNosel_QCD.push_back((TH1F*)ifile_QCD->Get(Form("Nosel_%s", cutNames[i].c_str())));
        hNosel_Gluino.push_back((TH1F*)ifile_Gluino->Get(Form("Nosel_%s", cutNames[i].c_str())));
        hNosel_MC.push_back((TH1F*)ifile_Wjets->Get(Form("Nosel_%s", cutNames[i].c_str())));
    }

    // Step 2: style + build MC sum
    for (size_t i = 0; i < cutNames.size(); i++) {
        cout << "Cut: " << cutNames[i] << endl;
        hNosel_MC[i]->Add(hNosel_TTbar[i]);
        hNosel_MC[i]->Add(hNosel_QCD[i]);

        // PUppiMET, Ptpseudo, TrkIso, EoverP, PtErr_over_PtPt, PtErr_over_Pt
        if (i==0 || i==1 || i==4 || i==5 || i==6 || i==7) {
            hNosel_TTbar[i]->Rebin(4);
            hNosel_TTbarSemiLep[i]->Rebin(4);
            hNosel_Wjets[i]->Rebin(4);
            hNosel_QCD[i]->Rebin(4);
            hNosel_Gluino[i]->Rebin(4);
            hNosel_JetMETdata[i]->Rebin(4);
        }
        // PFMiniIso
        if (i==3) {
            hNosel_TTbar[i]->Rebin(8);
            hNosel_TTbarSemiLep[i]->Rebin(8);
            hNosel_Wjets[i]->Rebin(8);
            hNosel_QCD[i]->Rebin(8);
            hNosel_Gluino[i]->Rebin(8);
            hNosel_JetMETdata[i]->Rebin(8);
        }

        hNosel_TTbar[i]->SetLineColor(kRed);
        hNosel_TTbar[i]->SetFillColorAlpha(kRed, 0.5);
        hNosel_TTbarSemiLep[i]->SetLineColor(kRed+2);
        hNosel_TTbarSemiLep[i]->SetFillColorAlpha(kRed+2, 0.5);
        hNosel_Wjets[i]->SetLineColor(kBlue-7);
        hNosel_Wjets[i]->SetFillColorAlpha(kBlue-7, 0.5);
        hNosel_QCD[i]->SetLineColor(kGreen);
        hNosel_QCD[i]->SetFillColorAlpha(kGreen, 0.5);
        hNosel_Gluino[i]->SetLineColor(kViolet-1);
        hNosel_Gluino[i]->SetMarkerStyle(22);
        hNosel_Gluino[i]->SetMarkerColor(kViolet-1);
        hNosel_JetMETdata[i]->SetLineColor(kBlack);
        hNosel_JetMETdata[i]->SetMarkerColor(kBlack);
        hNosel_JetMETdata[i]->SetMarkerStyle(20);
    }

    // Step 2 bis: if rescaled, scale MC to data
    if (isRescaled) {
        for (size_t i = 0; i < cutNames.size(); i++) {
            float intTot = hNosel_TTbar[i]->Integral() + hNosel_TTbarSemiLep[i]->Integral() + hNosel_Wjets[i]->Integral() + hNosel_QCD[i]->Integral();
            hNosel_TTbar[i]->Scale(hNosel_JetMETdata[i]->Integral()/intTot);
            hNosel_TTbarSemiLep[i]->Scale(hNosel_JetMETdata[i]->Integral()/intTot);
            hNosel_Wjets[i]->Scale(hNosel_JetMETdata[i]->Integral()/intTot);
            hNosel_QCD[i]->Scale(hNosel_JetMETdata[i]->Integral()/intTot);
            hNosel_Gluino[i]->Scale(hNosel_JetMETdata[i]->Integral()/hNosel_Gluino[i]->Integral());
        }
    }

    struct CutInfo {
        double threshold;
        bool keepRight;
        bool isSymmetric;
    };

    std::vector<CutInfo> cutInfos = {
        {150,    true,  false},  // PUppiMET      > 150
        {50,     true,  false},  // Pt_pseudo     > 50
        {2.4,    false, true },  // |eta|         < 2.4
        {0.02,   false, false},  // PFMiniIso     < 0.02
        {15,     false, false},  // TrkIso        < 15
        {0.3,    false, false},  // EoverP        < 0.3
        {0.0008, false, false},  // PtErr/PtPt    < 0.0008
        {1,      false, false},  // PtErr/Pt      < 1
        {0.3,    true,  false},  // Fpix          > 0.3
        {2.9784, true,  false},  // Ih
    };

    TLatex *tex = new TLatex(0.75, 0.91, "109 fb^{-1} (13.6 TeV)");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);

    TLatex *latex1 = new TLatex(0.16, 0.91, "#it{Private work (CMS simulation/data)}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);

    auto addOverflow = [](TH1F* h) {
        int n = h->GetNbinsX();
        h->SetBinContent(n, h->GetBinContent(n) + h->GetBinContent(n+1));
        h->SetBinError(n, std::sqrt(std::pow(h->GetBinError(n),2) + std::pow(h->GetBinError(n+1),2)));
        h->SetBinContent(n+1, 0);   // évite le double comptage
        h->SetBinError(n+1, 0);
    };

    // Step 3: Canvas for each hists
    std::vector<TCanvas*> canvases;
    for (size_t i = 0; i < cutNames.size(); i++) {

        TCanvas *c = new TCanvas(Form("c_Nosel_%s", cutNames[i].c_str()), Form("c_Nosel_%s", cutNames[i].c_str()), 800, 600);
        c->SetLeftMargin(0.16); c->SetBottomMargin(0.16);

        addOverflow(hNosel_QCD[i]);
        addOverflow(hNosel_TTbar[i]);
        addOverflow(hNosel_TTbarSemiLep[i]);
        addOverflow(hNosel_Wjets[i]);
        addOverflow(hNosel_JetMETdata[i]);
        addOverflow(hNosel_Gluino[i]);

        THStack *hs = new THStack(Form("hs_%s", cutNames[i].c_str()), "");
        hs->Add(hNosel_TTbar[i]);
        hs->Add(hNosel_TTbarSemiLep[i]);
        hs->Add(hNosel_QCD[i]);
        hs->Add(hNosel_Wjets[i]);

        double ymin = 0.001;
        double ymax = 1e12;
        double xmin = hNosel_TTbar[i]->GetBinLowEdge(1);
        double xmax = hNosel_TTbar[i]->GetBinLowEdge(hNosel_TTbar[i]->GetNbinsX()+1);

        if (i==9) xmax = 10;

        if (i==9) {
            hNosel_TTbar[i]->GetXaxis()->SetRangeUser(xmin, xmax);
            hNosel_TTbarSemiLep[i]->GetXaxis()->SetRangeUser(xmin, xmax);
            hNosel_Wjets[i]->GetXaxis()->SetRangeUser(xmin, xmax);
            hNosel_QCD[i]->GetXaxis()->SetRangeUser(xmin, xmax);
            hNosel_JetMETdata[i]->GetXaxis()->SetRangeUser(xmin, xmax);
            hNosel_Gluino[i]->GetXaxis()->SetRangeUser(xmin, xmax);
        }

        c->cd();
        TH1F *frame = (TH1F*)hNosel_TTbar[i]->Clone(Form("frame_%s", cutNames[i].c_str()));
        gStyle->SetOptStat(0);
        frame->Reset();
        frame->SetMinimum(ymin);
        frame->SetMaximum(ymax);
        frame->GetXaxis()->SetTitle(Xlabel[i].c_str());
        frame->GetYaxis()->SetTitle("Number of events");
        frame->GetYaxis()->SetTitleSize(0.06);
        frame->GetXaxis()->SetTitleSize(0.06);
        frame->GetXaxis()->SetTitleOffset(0.9);
        frame->GetYaxis()->SetTitleOffset(0.8);
        frame->GetXaxis()->SetLabelSize(0.05);
        frame->GetYaxis()->SetLabelSize(0.05);
        frame->GetXaxis()->SetRangeUser(xmin, xmax);

        frame->Draw("HIST");
        hs->Draw("HIST same");
        hNosel_JetMETdata[i]->Draw("E1 same");
        hNosel_Gluino[i]->Draw("E1 same");

        // ligne verticale à la valeur de la coupure
        double xcut = cutInfos[i].threshold;
        TLine *line = new TLine(xcut, ymin, xcut, ymax);
        line->SetLineColor(kBlack);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");

        // pour les coupures symétriques |var| < seuil, tracer aussi la ligne à -seuil
        if (cutInfos[i].isSymmetric) {
            TLine *lineNeg = new TLine(-xcut, ymin, -xcut, ymax);
            lineNeg->SetLineColor(kBlack);
            lineNeg->SetLineStyle(2);
            lineNeg->SetLineWidth(2);
            lineNeg->Draw("same");
        }

        tex->Draw();
        latex1->Draw();

        TLegend *legend = new TLegend(0.7, 0.6, 0.89, 0.89);
        if (i==2 || i==8) legend = new TLegend(0.67, 0.16, 0.89, 0.45);  // eta, Fpix
        legend->AddEntry(hNosel_TTbar[i],        "t#bar{t}#rightarrow2l2#nu",  "f");
        legend->AddEntry(hNosel_TTbarSemiLep[i], "t#bar{t}#rightarrowl#nu",    "f");
        legend->AddEntry(hNosel_Wjets[i],        "W(#rightarrow#mu#nu)+jets",  "f");
        legend->AddEntry(hNosel_QCD[i],          "QCD (#mu enriched)",         "f");
        legend->AddEntry(hNosel_JetMETdata[i],   "MET data",                   "lep");
        legend->AddEntry(hNosel_Gluino[i],       "#tilde{g} (m=2000 GeV)",     "lep");
        legend->SetBorderSize(0);

        TH1F *htemp = (TH1F*)hNosel_Wjets[i]->Clone("htemp");
        htemp->Add(hNosel_TTbar[i]);
        htemp->Add(hNosel_TTbarSemiLep[i]);
        htemp->Add(hNosel_QCD[i]);

        htemp->SetMinimum(ymin);
        htemp->SetMaximum(ymax);
        hNosel_JetMETdata[i]->SetMinimum(ymin);
        hNosel_JetMETdata[i]->SetMaximum(ymax);

        c->cd();
        gPad->RedrawAxis();
        c->SetLogy();
        c->SetTickx(1);
        c->SetTicky(1);
        hs->SetMinimum(ymin);
        hs->SetMaximum(ymax);
        c->Modified();
        legend->Draw();
        c->Update();

        TCanvas *cRatio = DrawWithCDF(hNosel_JetMETdata[i], htemp, hNosel_Gluino[i], c,
                                      Form("cCDF_%s", cutNames[i].c_str()), Xlabel[i], xmin, xmax);

        canvases.push_back(cRatio);
        delete htemp;
    }

    // Step 4: compute the efficiency for each cut
    auto computeEff = [](TH1F* h, double threshold, bool keepRight, bool isSymmetric = false) -> double {
        if (!h) return -1.0;
        int totalBins = h->GetNbinsX();
        double total  = h->Integral(1, totalBins);
        if (total == 0) return 0.0;
        if (isSymmetric) {
            int binLow  = h->FindBin(-threshold);
            int binHigh = h->FindBin(threshold);
            return h->Integral(binLow, binHigh) / total;
        }
        int cutBin = h->FindBin(threshold);
        double signal = keepRight ? h->Integral(cutBin, totalBins)
                                  : h->Integral(1, cutBin - 1);
        return signal / total;
    };

    // Print efficiencies
    std::cout << std::left
            << std::setw(25) << "Cut"
            << std::setw(12) << "TTbar"
            << std::setw(14) << "TTbarSemiLep"
            << std::setw(12) << "W+jets"
            << std::setw(12) << "QCD"
            << std::setw(12) << "ALL MC"
            << std::setw(12) << "Data"
            << std::endl;
    std::cout << std::string(60, '-') << std::endl;

    for (size_t i = 0; i < cutNames.size(); i++) {
        double eff_TTbar = computeEff(hNosel_TTbar[i],      cutInfos[i].threshold, cutInfos[i].keepRight, cutInfos[i].isSymmetric);
        double eff_TTbarSemiLep = computeEff(hNosel_TTbarSemiLep[i], cutInfos[i].threshold, cutInfos[i].keepRight, cutInfos[i].isSymmetric);
        double eff_Wjets = computeEff(hNosel_Wjets[i],      cutInfos[i].threshold, cutInfos[i].keepRight, cutInfos[i].isSymmetric);
        double eff_QCD   = computeEff(hNosel_QCD[i],        cutInfos[i].threshold, cutInfos[i].keepRight, cutInfos[i].isSymmetric);
        double eff_MC    = computeEff(hNosel_MC[i],         cutInfos[i].threshold, cutInfos[i].keepRight, cutInfos[i].isSymmetric);
        double eff_Data  = computeEff(hNosel_JetMETdata[i], cutInfos[i].threshold, cutInfos[i].keepRight, cutInfos[i].isSymmetric);

        std::cout << std::left  << std::setw(25) << cutNames[i]
                << std::fixed << std::setprecision(4)
                << std::setw(12) << eff_TTbar
                << std::setw(14) << eff_TTbarSemiLep
                << std::setw(12) << eff_Wjets
                << std::setw(12) << eff_QCD
                << std::setw(12) << eff_MC
                << std::setw(12) << eff_Data
                << std::endl;
    }

    // Step 5: save
    TString pdfName = Form("../outputDisplay/Nm1plots/NoselEff%s_ALL.pdf", isRescaled ? "_rescaled" : "");
    for (size_t i = 0; i < canvases.size(); i++) {
        if      (i == 0)                    canvases[i]->Print(pdfName + "(");
        else if (i == canvases.size() - 1)  canvases[i]->Print(pdfName + ")");
        else                                canvases[i]->Print(pdfName);
    }

    return;
}


void Nm1EffCutZeroed (bool isRescaled,
                      std::string TTbar,
                      std::string Wjets,
                      std::string QCD,
                      std::string JetMETdata,
                      std::string Gluino) {

    gErrorIgnoreLevel = kWarning;

    TFile *ifile_TTbar = new TFile(TTbar.c_str(), "READ");
    TFile *ifile_Wjets = new TFile(Wjets.c_str(), "READ");
    TFile *ifile_QCD = new TFile(QCD.c_str(), "READ");
    TFile *ifile_JetMETdata = new TFile(JetMETdata.c_str(), "READ");
    TFile *ifile_Gluino = new TFile(Gluino.c_str(), "READ");

    std::vector<TH1F*> hNM1_TTbar, hNM1_Wjets, hNM1_QCD, hNM1_MC, hNM1_JetMETdata, hNM1_Gluino;

    std::vector<string> cutNames = {"trigger", "METfilters", "PUppiMET", "Ptpseudo", "eta", "NOPH", "FOVH", "NOM", "HighPurity",
    "Chi2", "dZ", "dXY", "PFMiniIso", "TrkIso", "EoverP", "PtErr_over_PtPt", "Fpix", "PtErr_over_Pt", "Ih_StripOnly", "Ih_StripOnly_rescaled"};

    std::vector<string> Xlabel = {"or MET trigger", "MET filters", "PUppi MET", "p_{T} [GeV]", "#eta", "Nb pixel hits", "frac. valid hits", "Nb dE/dx", "HighPurity",
    "#chi^{2}/NDOF", "d_{z} [cm]", "d_{xy} [cm]", "I_{PF}^{rel}", "I_{dr03}^{trk}", "E/p", "#sigma_{p_{T}}/p_{T}^{2}", "F_{pixel}", "#sigma_{p_{T}}/p_{T}", "I_{h} [MeV/cm]", "I_{h} rescaled [MeV/cm]"};

    // Step 1: retrieve the Nm1 hists
    for (size_t i = 0; i < cutNames.size(); i++) {

        if (cutNames[i] == "Ih_StripOnly_rescaled") {
            hNM1_JetMETdata.push_back((TH1F*)ifile_JetMETdata->Get(Form("Nm1_event_%s", cutNames[i-1].c_str())));
            hNM1_TTbar.push_back((TH1F*)ifile_TTbar->Get(Form("Nm1_event_%s", cutNames[i].c_str())));
            hNM1_Wjets.push_back((TH1F*)ifile_Wjets->Get(Form("Nm1_event_%s", cutNames[i].c_str())));
            hNM1_QCD.push_back((TH1F*)ifile_QCD->Get(Form("Nm1_event_%s", cutNames[i].c_str())));
            hNM1_Gluino.push_back((TH1F*)ifile_Gluino->Get(Form("Nm1_event_%s", cutNames[i].c_str())));
            hNM1_MC.push_back((TH1F*)ifile_Wjets->Get(Form("Nm1_event_%s", cutNames[i].c_str())));
        }

        hNM1_JetMETdata.push_back((TH1F*)ifile_JetMETdata->Get(Form("Nm1_event_%s", cutNames[i].c_str())));
        hNM1_TTbar.push_back((TH1F*)ifile_TTbar->Get(Form("Nm1_event_%s", cutNames[i].c_str())));
        hNM1_Wjets.push_back((TH1F*)ifile_Wjets->Get(Form("Nm1_event_%s", cutNames[i].c_str())));
        hNM1_QCD.push_back((TH1F*)ifile_QCD->Get(Form("Nm1_event_%s", cutNames[i].c_str())));
        hNM1_Gluino.push_back((TH1F*)ifile_Gluino->Get(Form("Nm1_event_%s", cutNames[i].c_str())));
        hNM1_MC.push_back((TH1F*)ifile_Wjets->Get(Form("Nm1_event_%s", cutNames[i].c_str())));
    }

    // Step 2: style + build MC sum
    for (size_t i = 0; i < cutNames.size(); i++) {
        cout << "Cut: " << cutNames[i] << endl;
        hNM1_MC[i]->Add(hNM1_TTbar[i]);
        hNM1_MC[i]->Add(hNM1_QCD[i]);

        if (i==2 || i==3 || i==13 || i==14 || i==15 || i==17) {
            hNM1_TTbar[i]->Rebin(4);
            hNM1_Wjets[i]->Rebin(4);
            hNM1_QCD[i]->Rebin(4);
            hNM1_Gluino[i]->Rebin(4);
            hNM1_JetMETdata[i]->Rebin(4);
        }
        if (i==12) {
            hNM1_TTbar[i]->Rebin(8);
            hNM1_Wjets[i]->Rebin(8);
            hNM1_QCD[i]->Rebin(8);
            hNM1_Gluino[i]->Rebin(8);
            hNM1_JetMETdata[i]->Rebin(8);
        }


        hNM1_TTbar[i]->SetLineColor(kRed);
        hNM1_TTbar[i]->SetFillColorAlpha(kRed, 0.5);
        hNM1_Wjets[i]->SetLineColor(kBlue-7);
        hNM1_Wjets[i]->SetFillColorAlpha(kBlue-7, 0.5);
        hNM1_QCD[i]->SetLineColor(kGreen);
        hNM1_QCD[i]->SetFillColorAlpha(kGreen, 0.5);
        hNM1_Gluino[i]->SetLineColor(kViolet-1);
        hNM1_Gluino[i]->SetMarkerStyle(22);
        hNM1_Gluino[i]->SetMarkerColor(kViolet-1);
        hNM1_JetMETdata[i]->SetLineColor(kBlack);
        hNM1_JetMETdata[i]->SetMarkerColor(kBlack);
        hNM1_JetMETdata[i]->SetMarkerStyle(20);
    }

    // Step 2 bis: if rescaled, scale MC to data
    if (isRescaled) {
        for (size_t i = 0; i < cutNames.size(); i++) {
            float intTot = hNM1_TTbar[i]->Integral() + hNM1_Wjets[i]->Integral() + hNM1_QCD[i]->Integral();
            hNM1_TTbar[i]->Scale(hNM1_JetMETdata[i]->Integral()/intTot);
            hNM1_Wjets[i]->Scale(hNM1_JetMETdata[i]->Integral()/intTot);
            hNM1_QCD[i]->Scale(hNM1_JetMETdata[i]->Integral()/intTot);
            hNM1_Gluino[i]->Scale(hNM1_JetMETdata[i]->Integral()/hNM1_Gluino[i]->Integral());
        }
    }

    struct CutInfo {
        double threshold;
        bool keepRight;
        bool isSymmetric;
    };

    // Détail des coupures (cf. table cutflow) :
    // MET trigger (OR), MET filters, PUppiMET > 150 GeV, pT > 50 GeV, |eta| < 2.4,
    // N_pixel hits >= 2, f_valid hits > 0.8, N_dE/dx >= 10, HighPurity,
    // chi2/Ndof < 5, |dz| < 0.1 cm, |dxy| < 0.02 cm, I_PF^rel < 0.02, I_trk < 15 GeV,
    // E/p < 0.3, sigma_pT/pT^2 < 0.0008, F_pixel > 0.3, sigma_pT/pT < 1, Ih > C_mass
    std::vector<CutInfo> cutInfos = {
        {0.5,      true,  false},  // trigger
        {0.5,      true,  false},  // METfilters
        {150,    true,  false},  // PUppiMET > 150
        {50,     true,  false},  // Pt_pseudo     > 50
        {2.4,    false, true },  // |eta|         < 2.4
        {2,      true,  false},  // NOPH          >= 2
        {0.8,    true,  false},  // FOVH          > 0.8
        {10,     true,  false},  // NOM           >= 10
        {0.5,      true,  false},  // HighPurity
        {5.0,    false, false},  // Chi2          < 5
        {0.1,    false, true },  // |dZ|          < 0.1
        {0.02,   false, true },  // |dXY|         < 0.02
        {0.02,   false, false},  // PFMiniIso     < 0.02
        {15,     false, false},  // TrkIso        < 15
        {0.3,    false, false},  // EoverP        < 0.3
        {0.0008, false, false},  // PtErr/PtPt    < 0.0008
        {0.3,    true,  false},  // Fpix          > 0.3
        {1,      false, false},  // PtErr/Pt      < 1
        {2.9784, true,  false},  // Ih_StripOnly  > C_mass
        {2.9784, true,  false},  // Ih_StripOnly rescaled > C_mass
    };

    TLatex *tex = new TLatex(0.75, 0.91, "109 fb^{-1} (13.6 TeV)");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);

    TLatex *latex1 = new TLatex(0.16, 0.91, "#it{Private work (CMS simulation/data)}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);

    auto addOverflow = [](TH1F* h) {
        int n = h->GetNbinsX();
        h->SetBinContent(n, h->GetBinContent(n) + h->GetBinContent(n+1));
        h->SetBinError(n, std::sqrt(std::pow(h->GetBinError(n),2) + std::pow(h->GetBinError(n+1),2)));
        h->SetBinContent(n+1, 0);   // évite le double comptage
        h->SetBinError(n+1, 0);
    };

    // Met à zéro tous les bins en dehors de la région acceptée par la coupure.
    // Conventions identiques à computeEff : keepRight garde [FindBin(thr), n],
    // keepLeft garde [1, FindBin(thr)-1], symétrique garde [FindBin(-thr), FindBin(thr)].
    auto zeroOutsideCut = [](TH1F* h, double threshold, bool keepRight, bool isSymmetric) {
        if (!h) return;
        int n = h->GetNbinsX();
        if (isSymmetric) {
            int binLow  = h->FindBin(-threshold);
            int binHigh = h->FindBin(threshold);
            for (int b = 0; b <= n+1; b++) {
                if (b < binLow || b > binHigh) {
                    h->SetBinContent(b, 0);
                    h->SetBinError(b, 0);
                }
            }
            return;
        }
        int cutBin = h->FindBin(threshold);
        for (int b = 0; b <= n+1; b++) {
            bool keep = keepRight ? (b >= cutBin) : (b <= cutBin - 1);
            if (!keep) {
                h->SetBinContent(b, 0);
                h->SetBinError(b, 0);
            }
        }
    };

    // Step 3 (déplacé avant les canvas) : compute the efficiency for each cut
    // NB : doit être fait AVANT la mise à zéro des bins hors coupure,
    // sinon toutes les efficacités valent trivialement 1.
    auto computeEff = [](TH1F* h, double threshold, bool keepRight, bool isSymmetric = false) -> double {
        if (!h) return -1.0;
        int totalBins = h->GetNbinsX();
        double total  = h->Integral(1, totalBins);
        if (total == 0) return 0.0;
        if (isSymmetric) {
            int binLow  = h->FindBin(-threshold);
            int binHigh = h->FindBin(threshold);
            return h->Integral(binLow, binHigh) / total;
        }
        int cutBin = h->FindBin(threshold);
        double signal = keepRight ? h->Integral(cutBin, totalBins)
                                  : h->Integral(1, cutBin - 1);
        return signal / total;
    };

    // Print efficiencies
    std::cout << std::left
            << std::setw(25) << "Cut"
            << std::setw(12) << "TTbar"
            << std::setw(12) << "W+jets"
            << std::setw(12) << "QCD"
            << std::setw(12) << "ALL MC"
            << std::setw(12) << "Data"
            << std::endl;
    std::cout << std::string(60, '-') << std::endl;

    for (size_t i = 0; i < cutNames.size(); i++) {
        double eff_TTbar = computeEff(hNM1_TTbar[i],      cutInfos[i].threshold, cutInfos[i].keepRight, cutInfos[i].isSymmetric);
        double eff_Wjets = computeEff(hNM1_Wjets[i],      cutInfos[i].threshold, cutInfos[i].keepRight, cutInfos[i].isSymmetric);
        double eff_QCD   = computeEff(hNM1_QCD[i],        cutInfos[i].threshold, cutInfos[i].keepRight, cutInfos[i].isSymmetric);
        double eff_MC    = computeEff(hNM1_MC[i],         cutInfos[i].threshold, cutInfos[i].keepRight, cutInfos[i].isSymmetric);
        double eff_Data  = computeEff(hNM1_JetMETdata[i], cutInfos[i].threshold, cutInfos[i].keepRight, cutInfos[i].isSymmetric);

        std::cout << std::left  << std::setw(25) << cutNames[i]
                << std::fixed << std::setprecision(4)
                << std::setw(12) << eff_TTbar
                << std::setw(12) << eff_Wjets
                << std::setw(12) << eff_QCD
                << std::setw(12) << eff_MC
                << std::setw(12) << eff_Data
                << std::endl;
    }

    // Step 4: Canvas for each hists (bins hors coupure mis à zéro)
    std::vector<TCanvas*> canvases;
    for (size_t i = 0; i < cutNames.size(); i++) {

        TCanvas *c = new TCanvas(Form("c_Nm1_%s", cutNames[i].c_str()), Form("c_Nm1_%s", cutNames[i].c_str()), 800, 600);
        c->SetLeftMargin(0.16); c->SetBottomMargin(0.16);

        addOverflow(hNM1_QCD[i]);
        addOverflow(hNM1_TTbar[i]);
        addOverflow(hNM1_Wjets[i]);
        addOverflow(hNM1_JetMETdata[i]);
        addOverflow(hNM1_Gluino[i]);

        // ----- subtilité : on met à zéro les bins qui ne passent pas la coupure -----
        zeroOutsideCut(hNM1_QCD[i],        cutInfos[i].threshold, cutInfos[i].keepRight, cutInfos[i].isSymmetric);
        zeroOutsideCut(hNM1_TTbar[i],      cutInfos[i].threshold, cutInfos[i].keepRight, cutInfos[i].isSymmetric);
        zeroOutsideCut(hNM1_Wjets[i],      cutInfos[i].threshold, cutInfos[i].keepRight, cutInfos[i].isSymmetric);
        zeroOutsideCut(hNM1_JetMETdata[i], cutInfos[i].threshold, cutInfos[i].keepRight, cutInfos[i].isSymmetric);
        zeroOutsideCut(hNM1_Gluino[i],     cutInfos[i].threshold, cutInfos[i].keepRight, cutInfos[i].isSymmetric);
        // -----------------------------------------------------------------------------

        THStack *hs = new THStack(Form("hs_%s", cutNames[i].c_str()), "");
        hs->Add(hNM1_QCD[i]);
        hs->Add(hNM1_TTbar[i]);
        hs->Add(hNM1_Wjets[i]);

        double ymin = 0.01;
        double ymax = 6e6;
        double xmin = hNM1_TTbar[i]->GetBinLowEdge(1);
        double xmax = hNM1_TTbar[i]->GetBinLowEdge(hNM1_TTbar[i]->GetNbinsX()+1);

        if (i==6) xmin = 0.75;

        c->cd();
        TH1F *frame = (TH1F*)hNM1_TTbar[i]->Clone(Form("frame_%s", cutNames[i].c_str()));
        gStyle->SetOptStat(0);
        frame->Reset();
        frame->SetMinimum(ymin);
        frame->SetMaximum(ymax);
        frame->GetXaxis()->SetTitle(Xlabel[i].c_str());
        frame->GetYaxis()->SetTitle("Number of events");
        frame->GetYaxis()->SetTitleSize(0.06);
        frame->GetXaxis()->SetTitleSize(0.06);
        frame->GetXaxis()->SetTitleOffset(0.9);
        frame->GetYaxis()->SetTitleOffset(0.8);
        frame->GetXaxis()->SetLabelSize(0.05);
        frame->GetYaxis()->SetLabelSize(0.05);
        frame->GetXaxis()->SetRangeUser(xmin, xmax);

        frame->Draw("HIST");
        hs->Draw("HIST same");
        hNM1_JetMETdata[i]->Draw("E1 same");
        hNM1_Gluino[i]->Draw("E1 same");


        // ligne verticale à la valeur de la coupure
        double xcut = cutInfos[i].threshold;
        TLine *line = new TLine(xcut, ymin, xcut, ymax);
        line->SetLineColor(kBlack);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");

        // pour les coupures symétriques |var| < seuil, tracer aussi la ligne à -seuil
        if (cutInfos[i].isSymmetric) {
            TLine *lineNeg = new TLine(-xcut, ymin, -xcut, ymax);
            lineNeg->SetLineColor(kBlack);
            lineNeg->SetLineStyle(2);
            lineNeg->SetLineWidth(2);
            lineNeg->Draw("same");
        }

        tex->Draw();
        latex1->Draw();


        TLegend *legend = new TLegend(0.7, 0.6, 0.89, 0.89);
        if (i==0 || i==1 || i==6 || i==4 || i==8 || i==16) legend = new TLegend(0.67, 0.16, 0.89, 0.45);
        legend->AddEntry(hNM1_TTbar[i],     "TTbar",     "f");
        legend->AddEntry(hNM1_Wjets[i],     "W(#rightarrow#mu#nu)+jets",    "f");
        legend->AddEntry(hNM1_QCD[i],       "QCD (#mu enriched)",       "f");
        legend->AddEntry(hNM1_JetMETdata[i],"MET data","lep");
        legend->AddEntry(hNM1_Gluino[i],     "#tilde{g} (m=2000 GeV)", "lep");
        legend->SetBorderSize(0);

        TH1F *htemp = (TH1F*)hNM1_Wjets[i]->Clone("htemp");
        htemp->Add(hNM1_TTbar[i]);
        htemp->Add(hNM1_QCD[i]);

        htemp->SetMinimum(ymin);
        htemp->SetMaximum(ymax);
        hNM1_JetMETdata[i]->SetMinimum(ymin);
        hNM1_JetMETdata[i]->SetMaximum(ymax);

        c->cd();
        gPad->RedrawAxis();
        c->SetLogy();
        c->SetTickx(1);
        c->SetTicky(1);
        hs->SetMinimum(ymin);
        hs->SetMaximum(ymax);
        c->Modified();
        legend->Draw();
        c->Update();

        TCanvas *cRatio = DrawWithCDF(hNM1_JetMETdata[i], htemp, hNM1_Gluino[i], c,
                                      Form("cCDF_%s", cutNames[i].c_str()), Xlabel[i], xmin, xmax);

        canvases.push_back(cRatio);
        delete htemp;
    }

    // Step 5: save
    TString pdfName = Form("../outputDisplay/Nm1plots/Nm1EffCutZeroed%s_ALL.pdf", isRescaled ? "_rescaled" : "");
    for (size_t i = 0; i < canvases.size(); i++) {
        if      (i == 0)                    canvases[i]->Print(pdfName + "(");
        else if (i == canvases.size() - 1)  canvases[i]->Print(pdfName + ")");
        else                                canvases[i]->Print(pdfName);
    }

    return;
}

void PlotPthatQCD(std::string sampleW  = "../output/QCD2024_V16/QCD2024_mu_V16p2_weighted.root",
                  std::string sampleNW = "../output/QCD2024_V16/QCD2024_mu_V16p2.root") {

    gErrorIgnoreLevel = kWarning;
    gStyle->SetOptStat(0);

    TFile *ifileW  = new TFile(sampleW.c_str(),  "READ");
    TFile *ifileNW = new TFile(sampleNW.c_str(), "READ");
    if (!ifileW || ifileW->IsZombie()) {
        std::cerr << "Cannot open file: " << sampleW << std::endl;
        return;
    }
    if (!ifileNW || ifileNW->IsZombie()) {
        std::cerr << "Cannot open file: " << sampleNW << std::endl;
        return;
    }

    TH1F *hW  = (TH1F*)ifileW->Get("Nosel_PthatQCD");
    TH1F *hNW = (TH1F*)ifileNW->Get("Nosel_PthatQCD");
    if (!hW) {
        std::cerr << "Histogram Nosel_PthatQCD not found in " << sampleW << std::endl;
        return;
    }
    if (!hNW) {
        std::cerr << "Histogram Nosel_PthatQCD not found in " << sampleNW << std::endl;
        return;
    }

    TCanvas *c = new TCanvas("c_PthatQCD", "c_PthatQCD", 800, 600);
    c->SetLeftMargin(0.16);
    c->SetBottomMargin(0.16);
    c->SetLogy();
    c->SetTickx(1);
    c->SetTicky(1);

    hNW->SetLineColor(kBlack);
    hNW->SetLineWidth(2);

    hW->SetLineColor(kGreen+2);
    hW->SetLineWidth(2);
    hW->GetXaxis()->SetTitle("#hat{p}_{T} [GeV]");
    hW->GetYaxis()->SetTitle("Number of events");
    hW->GetXaxis()->SetTitleSize(0.06);
    hW->GetYaxis()->SetTitleSize(0.06);
    hW->GetXaxis()->SetTitleOffset(0.9);
    hW->GetYaxis()->SetTitleOffset(0.8);
    hW->GetXaxis()->SetLabelSize(0.05);
    hW->GetYaxis()->SetLabelSize(0.05);

    double ymax = std::max(hW->GetMaximum(), hNW->GetMaximum()) * 10;
    hW->SetMaximum(ymax);
    hW->SetMinimum(0.1);

    hW->Draw("HIST");
    hNW->Draw("HIST same");

    TLatex *tex = new TLatex(0.68, 0.91, "109 fb^{-1} (13.6 TeV)");
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);
    tex->Draw();

    TLatex *latex1 = new TLatex(0.16, 0.91, "#it{Private work (CMS simulation)}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);
    latex1->Draw();

    TLegend *legend = new TLegend(0.6, 0.7, 0.89, 0.89);
    legend->AddEntry(hW,  "QCD weighted",     "l");
    legend->AddEntry(hNW, "QCD non-weighted", "l");
    legend->SetBorderSize(0);
    legend->Draw();

    gPad->RedrawAxis();
    c->Update();

    c->Print("../outputDisplay/Nosel_PthatQCD.pdf");

    return;
}

void TriggerEfficiency_VsBeta(const char *ofilename = "TriggEffBeta") {

    gErrorIgnoreLevel = kError;
    gStyle->SetOptStat(0);

    auto fileName = [&](int m) {
        return Form("../output/Gluino_V19/Gluino_Run3_MET_madgraph_%d_V19p10.root", m);
    };

    std::vector<int> masses = {2000, 2400, 2600};

    const std::string hDenomName = "Nosel_GenHSCP_beta";
    const std::string hNumName   = "Nosel_GenHSCP_beta_ifHLTMu50";

    const int rebin = 1;

    auto colorFor = [](int j) {
        if (j == 0) return (int)(kOrange+8);
        if (j == 1) return (int)(kViolet+1);
        if (j == 2) return (int)(kGreen-3);
        return 1;
    };

    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);

    std::vector<TEfficiency*> hEff;
    std::vector<TH1F*>        hBeta;
    std::vector<TH1*>         hCDF;
    std::vector<int>          hMassVal;

    for (int m : masses) {

        TFile *f = new TFile(fileName(m), "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "Warning: cannot open " << fileName(m) << ", skipping mass " << m << std::endl;
            continue;
        }

        TH1 *hDenom = (TH1*)f->Get(hDenomName.c_str());
        TH1 *hNum   = (TH1*)f->Get(hNumName.c_str());
        if (!hDenom || !hNum) {
            std::cerr << "Warning: missing histos in " << fileName(m) << ", skipping mass " << m << std::endl;
            f->Close();
            continue;
        }

        TH1F *hD = (TH1F*)hDenom->Clone(Form("hD_m%d", m));
        TH1F *hN = (TH1F*)hNum->Clone(Form("hN_m%d", m));
        hD->SetDirectory(0);
        hN->SetDirectory(0);
        if (rebin > 1) { hD->Rebin(rebin); hN->Rebin(rebin); }

        // copie du denominateur AVANT le garde-fou, pour la distribution brute
        TH1F *hDist = (TH1F*)hD->Clone(Form("hBeta_m%d", m));
        hDist->SetDirectory(0);
        if (hDist->Integral() > 0) hDist->Scale(1. / hDist->Integral());

        // CDF construite sur la distribution deja normalisee -> va de 0 a 1
        TH1 *hc = hDist->GetCumulative();
        hc->SetName(Form("hBetaCDF_m%d", m));
        hc->SetDirectory(0);

        // securite : num <= denom bin par bin (poids PU peuvent casser ca)
        for (int b = 1; b <= hD->GetNbinsX(); ++b) {
            if (hN->GetBinContent(b) > hD->GetBinContent(b))
                hN->SetBinContent(b, hD->GetBinContent(b));
        }

        if (!TEfficiency::CheckConsistency(*hN, *hD)) {
            std::cerr << "Warning: inconsistent histos for mass " << m << ", skipping" << std::endl;
            f->Close();
            continue;
        }

        TEfficiency *te = new TEfficiency(*hN, *hD);
        te->SetName(Form("teBeta_m%d", m));
        te->SetStatisticOption(TEfficiency::kFCP);

        f->Close();

        hEff.push_back(te);
        hBeta.push_back(hDist);
        hCDF.push_back(hc);
        hMassVal.push_back(m);
    }

    if (hEff.empty()) return;

    int nCurves = (int)hEff.size();

    double distMax = 0.;
    for (int j = 0; j < nCurves; ++j) {
        int ci = colorFor(j);

        hEff[j]->SetMarkerStyle(20+j);
        hEff[j]->SetMarkerColor(ci);
        hEff[j]->SetLineColor(ci);
        hEff[j]->SetLineWidth(2);
        hEff[j]->SetMarkerSize(1.1);

        hBeta[j]->SetLineColor(ci);
        hBeta[j]->SetLineWidth(2);
        hBeta[j]->SetMarkerColor(ci);
        hBeta[j]->SetMarkerStyle(20+j);
        hBeta[j]->SetFillStyle(0);

        hCDF[j]->SetLineColor(ci);
        hCDF[j]->SetLineWidth(2);
        hCDF[j]->SetMarkerColor(ci);
        hCDF[j]->SetMarkerStyle(20+j);
        hCDF[j]->SetMarkerSize(0.7);
        hCDF[j]->SetFillStyle(0);

        distMax = std::max(distMax, hBeta[j]->GetMaximum());
    }

    auto styleFrame = [&](TH1F *fr, const char *ytitle, double ymin, double ymax) {
        fr->SetDirectory(0);
        fr->GetXaxis()->SetTitle("#beta_{gen}");
        fr->GetYaxis()->SetTitle(ytitle);
        fr->GetYaxis()->SetTitleSize(0.055);
        fr->GetXaxis()->SetTitleSize(0.06);
        fr->GetYaxis()->SetTitleOffset(1.2);
        fr->GetXaxis()->SetTitleOffset(1.2);
        fr->GetYaxis()->SetLabelSize(0.05);
        fr->GetXaxis()->SetLabelSize(0.045);
        fr->SetMinimum(ymin);
        fr->SetMaximum(ymax);
    };

    auto makeLeg = [&]() {
        TLegend *l = new TLegend(0.18, 0.55, 0.45, 0.89);
        l->SetBorderSize(0);
        for (int j = 0; j < nCurves; ++j)
            l->AddEntry(hEff[j], Form("m_{#tilde{g}} = %d GeV", hMassVal[j]), "lep");
        return l;
    };

    // ================= Canvas 1 : distribution de beta + CDF =================
    TCanvas *c1 = new TCanvas("c_BetaDistribution", "c_BetaDistribution", 800, 600);

    TPad *pad1 = new TPad("padBeta1", "padBeta1", 0.0, 0.30, 1.0, 1.0);
    pad1->SetLeftMargin(0.16);
    pad1->SetRightMargin(0.04);
    pad1->SetBottomMargin(0.20);
    pad1->SetGrid();
    pad1->Draw();

    TPad *pad2 = new TPad("padBeta2", "padBeta2", 0.0, 0.0, 1.0, 0.315);
    pad2->SetLeftMargin(0.16);
    pad2->SetRightMargin(0.04);
    pad2->SetBottomMargin(0.35);
    pad2->SetTopMargin(0.04);
    pad2->Draw();

    // --- pad du haut : distributions ---
    pad1->cd();
    TH1F *frameDist = new TH1F("frameBetaDist", "", 100, 0., 1.);
    styleFrame(frameDist, "Normalised entries", 0.0, 1.35 * distMax);
    frameDist->GetXaxis()->SetTitle("#beta_{gen}");
    frameDist->GetXaxis()->SetLabelSize(0.065);   // axe X porte par le pad du bas
    frameDist->GetXaxis()->SetTitleSize(0.075);
    frameDist->GetYaxis()->SetTitleSize(0.075);
    frameDist->GetYaxis()->SetTitleOffset(0.9);
    frameDist->GetXaxis()->SetTitleOffset(0.9);
    frameDist->GetYaxis()->SetLabelSize(0.065);
    frameDist->Draw();
    for (int j = 0; j < nCurves; ++j) hBeta[j]->Draw("HIST same");
    TLegend *leg1 = makeLeg();
    leg1->Draw();
    latex1->Draw();

    // --- pad du bas : CDF ---
    pad2->cd();
    TH1F *frameCDF = new TH1F("frameBetaCDF", "", 100, 0., 1.);
    frameCDF->SetDirectory(0);
    frameCDF->SetTitle("");
    frameCDF->GetYaxis()->SetTitle("CDF");
    frameCDF->GetXaxis()->SetTitle("#beta_{gen}");
    frameCDF->GetYaxis()->SetRangeUser(0., 1.);
    frameCDF->GetYaxis()->SetNdivisions(505);
    frameCDF->GetYaxis()->SetTitleFont(43);
    frameCDF->GetXaxis()->SetTitleFont(43);
    frameCDF->GetYaxis()->SetLabelFont(43);
    frameCDF->GetXaxis()->SetLabelFont(43);
    frameCDF->GetYaxis()->SetTitleSize(30);
    frameCDF->GetXaxis()->SetTitleSize(28);
    frameCDF->GetYaxis()->SetLabelSize(28);
    frameCDF->GetXaxis()->SetLabelSize(28);
    frameCDF->GetYaxis()->SetTitleOffset(1.3);
    frameCDF->GetXaxis()->SetTitleOffset(0.9);
    frameCDF->Draw();

    for (int j = 0; j < nCurves; ++j) hCDF[j]->Draw("hist P same");

    TLine *line50 = new TLine(0., 0.5, 1., 0.5);
    line50->SetLineColor(kBlack);
    line50->SetLineStyle(2);
    line50->Draw("same");

    c1->Modified(); c1->Update();
    c1->SaveAs(Form("../outputDisplay/TriggEff/c_BetaDistribution__%s.pdf", ofilename));

    // ================= Canvas 2 : efficacite =================
    TCanvas *c2 = new TCanvas("c_TriggerEfficiencyBeta", "c_TriggerEfficiencyBeta", 800, 600);
    c2->cd();
    c2->SetGrid();
    c2->SetLeftMargin(0.16);
    c2->SetRightMargin(0.04);
    c2->SetBottomMargin(0.19);

    TH1F *frameEff = new TH1F("frameBetaEff", "", 100, 0., 1.);
    styleFrame(frameEff, "Trigger efficiency (HLT_Mu50)", 0.0, 0.4);
    frameEff->Draw();
    for (int j = 0; j < nCurves; ++j) hEff[j]->Draw("E1 same");
    TLegend *leg2 = makeLeg();
    leg2->Draw();
    latex1->Draw();

    c2->Modified(); c2->Update();
    c2->SaveAs(Form("../outputDisplay/TriggEff/c_TriggerEfficiencyBeta__%s.pdf", ofilename));

    // ================= versions "Private work" =================
    latex1->SetTitle("#it{Private work (CMS simulation)}");

    c1->Modified(); c1->Update();
    c1->SaveAs(Form("../outputDisplay/TriggEff/c_BetaDistribution__%s_bis.pdf", ofilename));

    c2->Modified(); c2->Update();
    c2->SaveAs(Form("../outputDisplay/TriggEff/c_TriggerEfficiencyBeta__%s_bis.pdf", ofilename));

    return;
}


void Ihand1oP_fits_bkg() {

    // ---------------------------------------------------------------------
    // Configuration
    // ---------------------------------------------------------------------
    const TString inputPath =
        "/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/massSpectrum_bckgPrediction/DebugFit/Fits_HistForBkg_MC_8fp9_OldFit_Eta2p4_161_2.978400.root";

    const TString ihHistName = "proj_ih_eta10";
    const TString opHistName = "forfit_p_eta16";

    const double ihFitMax  = 5.5;    // upper edge of the Gaussian fit range
    const double ihDrawMax = 5.5;    // x-axis display range for Ih
    const double opDrawMax = 90.;  // x-axis display range for 1/p
    const double opDrawFitMax = 40;// range over which the 1/p fit is drawn

    const TString outDir = "../outputDisplay";

    // ---------------------------------------------------------------------
    // Input
    // ---------------------------------------------------------------------
    TFile *ifile = TFile::Open(inputPath, "READ");
    if (!ifile || ifile->IsZombie()) {
        Error("Ihand1oP_fits_bkg", "cannot open %s", inputPath.Data());
        return;
    }

    TH1 *Ih_new = (TH1*)ifile->Get(ihHistName);
    TH1 *oP_new = (TH1*)ifile->Get(opHistName);

    if (!Ih_new || !oP_new) {
        Error("Ihand1oP_fits_bkg", "missing histogram(s): %s %s",
              Ih_new ? "" : ihHistName.Data(), oP_new ? "" : opHistName.Data());
        ifile->ls();
        return;
    }
    Ih_new->SetDirectory(0);
    oP_new->SetDirectory(0);

    // ---------------------------------------------------------------------
    // 1/p fit
    // ---------------------------------------------------------------------
    float rangemax_p = 0.9 * oP_new->GetBinCenter(oP_new->GetMaximumBin());

    TF1 f_p_new("f_p_new", "[0]*([1]+erf((log(x)-[2])/[3]))", 0, rangemax_p);
    f_p_new.SetParameter(0, 1);
    f_p_new.FixParameter(1, 1.0);
    f_p_new.SetParameter(2, 3.50116e+00);
    f_p_new.SetParameter(3, 0.60152e+00);

    oP_new->Fit(&f_p_new, "RME", "", 0, rangemax_p);

    // ---------------------------------------------------------------------
    // Ih fit
    // ---------------------------------------------------------------------
    float max_ih = Ih_new->GetBinCenter(Ih_new->GetMaximumBin());

    float start_fit = 1.1 * max_ih;
    int lastBinContent = Ih_new->GetNbinsX();
    while (lastBinContent > 1 && Ih_new->GetBinContent(lastBinContent) == 0) lastBinContent--;
    if (start_fit > Ih_new->GetBinCenter(lastBinContent)) start_fit = max_ih;

    TF1 f_ih_new("f_ih_new", "gaus", start_fit, ihFitMax);
    f_ih_new.SetParameter(0, 0.5 * Ih_new->Integral());
    f_ih_new.SetParameter(1, max_ih);
    f_ih_new.SetParameter(2, Ih_new->GetStdDev());

    Ih_new->Fit(&f_ih_new, "RL", "", start_fit, ihFitMax);

    // ---------------------------------------------------------------------
    // Styling
    // ---------------------------------------------------------------------
    Ih_new->GetYaxis()->SetTitleSize(0.06);
    Ih_new->GetXaxis()->SetTitleSize(0.06);
    Ih_new->GetXaxis()->SetTitleOffset(0.9);
    Ih_new->GetYaxis()->SetTitleOffset(1);
    Ih_new->GetXaxis()->SetLabelSize(0.05);
    Ih_new->GetYaxis()->SetLabelSize(0.05);
    Ih_new->GetXaxis()->SetTitle("I_{h} [MeV/cm]");
    Ih_new->GetYaxis()->SetTitle("Number of tracks");
    Ih_new->SetLineColor(kBlack);
    Ih_new->SetMarkerStyle(20);
    Ih_new->SetMarkerColor(kBlack);
    Ih_new->GetXaxis()->SetRangeUser(0, ihDrawMax);

    oP_new->GetYaxis()->SetTitleSize(0.06);
    oP_new->GetXaxis()->SetTitleSize(0.06);
    oP_new->GetXaxis()->SetTitleOffset(0.9);
    oP_new->GetYaxis()->SetTitleOffset(1);
    oP_new->GetXaxis()->SetLabelSize(0.05);
    oP_new->GetYaxis()->SetLabelSize(0.05);
    oP_new->GetXaxis()->SetTitle("10^{4}/p [GeV^{-1}]");
    oP_new->GetYaxis()->SetTitle("Number of tracks (normalised)");
    oP_new->SetLineColor(kBlack);
    oP_new->SetMarkerStyle(20);
    oP_new->SetMarkerColor(kBlack);
    oP_new->GetXaxis()->SetRangeUser(0, opDrawMax);

    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);

    TLatex *etatex = new TLatex(0.71, 0.91, "#scale[1.3]{#bf{1.20#leq#eta<1.44}}");
    etatex->SetNDC();
    etatex->SetTextFont(42);
    etatex->SetTextSize(0.04);

    

    // ---------------------------------------------------------------------
    // Ih canvas
    // ---------------------------------------------------------------------
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    c1->cd();
    c1->SetLeftMargin(0.16); c1->SetBottomMargin(0.16);
    gStyle->SetOptStat(0);
    Ih_new->Draw("EO");
    latex1->Draw();
    etatex->Draw();
    c1->SetLogy();
    gStyle->SetOptFit(1);
    f_ih_new.SetLineColor(kRed);
    f_ih_new.SetLineWidth(2);
    f_ih_new.Draw("same");

    c1->Update();
    TPaveStats *st = (TPaveStats*)Ih_new->FindObject("stats");
    if (st) {
        st->SetX1NDC(0.47); // new x start position
        st->SetX2NDC(0.89); // new x end position
        st->SetY1NDC(0.6); // new y start position
        st->SetY2NDC(0.89); // new y end position
        st->SetTextSize(0.035);
    }
    c1->SaveAs(outDir + "/Ih_bkgfit_" + ihHistName + ".pdf");
    latex1->SetTitle("#it{Private work (CMS data)}");
    c1->Modified();
    c1->Update();
    c1->SaveAs(outDir + "/Ih_bkgfit_" + ihHistName + "_bis.pdf");
    latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Work in progress}");

    // ---------------------------------------------------------------------
    // 1/p canvas
    // ---------------------------------------------------------------------
    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    c2->cd();
    c2->SetLeftMargin(0.16); c2->SetBottomMargin(0.16);
    gStyle->SetOptStat(0);
    oP_new->Draw("E0");
    latex1->Draw();
    etatex->Draw();
    c2->SetLogy();
    gStyle->SetOptFit(1);
    f_p_new.SetLineColor(kRed);
    f_p_new.SetLineWidth(2);
    f_p_new.SetRange(0, opDrawFitMax);
    f_p_new.Draw("same");

    c2->Update();
    TPaveStats *st2 = (TPaveStats*)oP_new->FindObject("stats");
    if (st2) {
        st2->SetX1NDC(0.3); // new x start position
        st2->SetX2NDC(0.6); // new x end position
        st2->SetY1NDC(0.2); // new y start position
        st2->SetY2NDC(0.5); // new y end position
        st2->SetTextSize(0.035);
    }
    c2->SaveAs(outDir + "/oP_bkgfit_" + opHistName + ".pdf");
    latex1->SetTitle("#it{Private work (CMS data)}");
    c2->Modified();
    c2->Update();
    c2->SaveAs(outDir + "/oP_bkgfit_" + opHistName + "_bis.pdf");

    return;
}

// ---------------------------------------------------------------------------
// Overlay of background distributions in three eta slices, with a cumulative
// distribution panel underneath. Fit functions stored in the input file are
// stripped before drawing.
//
//   root -l etaSlices_bkg.C
// ---------------------------------------------------------------------------

const int nSlices = 3;
const TString etaLabel[nSlices] = { "0.00 #leq #eta < 0.24",
                                    "0.96 #leq #eta < 1.20",
                                    "2.16 #leq #eta < 2.40" };
const int sliceColor[nSlices]  = { kBlack, kRed, kBlue+1 };
const int sliceMarker[nSlices] = { 20, 21, 22 };
const TString outDir = "../outputDisplay";

bool loadEtaSlices (TFile *ifile,
                   const TString histName[nSlices],
                   TH1 *h[nSlices],
                   const TString tag,
                   bool normalise = true) {
    for (int s = 0; s < nSlices; ++s) {
        TH1 *hIn = (TH1*)ifile->Get(histName[s]);
        if (!hIn) {
            Error("loadEtaSlices", "missing histogram: %s", histName[s].Data());
            ifile->ls();
            return false;
        }
        h[s] = (TH1*)hIn->Clone(Form("%s_%s_clone", histName[s].Data(), tag.Data()));
        h[s]->SetDirectory(0);

        // drop any TF1 / TPaveStats attached to the histogram in the file
        if (h[s]->GetListOfFunctions()) h[s]->GetListOfFunctions()->Delete();
        h[s]->SetStats(0);

        if (normalise && h[s]->Integral() > 0) h[s]->Scale(1./h[s]->Integral());

        h[s]->SetLineColor(sliceColor[s]);
        h[s]->SetMarkerColor(sliceColor[s]);
        h[s]->SetMarkerStyle(sliceMarker[s]);
        h[s]->SetLineWidth(2);
    }
    return true;
}


TCanvas *drawEtaSlices (TFile *ifile, const TString histName[nSlices], const TString xTitle, double xDrawMin, double xDrawMax, double yHeadroom,
                       const TString canvasName, const TString outName, bool normalise = true, bool islogY = true) {

    TH1 *h[nSlices] = {nullptr};
    if (!loadEtaSlices(ifile, histName, h, canvasName, normalise)) return nullptr;

    h[0]->SetTitle("");
    h[0]->GetYaxis()->SetTitleSize(0.06);
    h[0]->GetYaxis()->SetTitleOffset(0.9);
    h[0]->GetYaxis()->SetLabelSize(0.05);
    h[0]->GetYaxis()->SetTitle(normalise ? "Number of tracks (normalised)" : "Number of tracks");
    h[0]->GetXaxis()->SetTitle(xTitle);
    h[0]->GetXaxis()->SetTitleSize(0.06);
    h[0]->GetXaxis()->SetTitleOffset(1.0);
    h[0]->GetXaxis()->SetLabelSize(0.05);
    h[0]->GetXaxis()->SetRangeUser(xDrawMin, xDrawMax);

    
    double yMax = 0., yMin = 1.e30;
    for (int s = 0; s < nSlices; ++s) {
        const int b1 = h[s]->FindBin(xDrawMin);
        const int b2 = h[s]->FindBin(xDrawMax);
        for (int b = b1; b <= b2; ++b) {
            const double y = h[s]->GetBinContent(b);
            if (y > yMax)          yMax = y;
            if (y > 0 && y < yMin) yMin = y;
        }
    }
    h[0]->SetMaximum(yHeadroom * yMax);
    h[0]->SetMinimum(islogY ? 0.5 * yMin : 0.);

    TH1 *hCDF[nSlices] = {nullptr};
    for (int s = 0; s < nSlices; ++s) {
        TH1 *hc = (TH1*)h[s]->Clone(Form("%s_forcdf", h[s]->GetName()));
        if (!normalise && hc->Integral() > 0) hc->Scale(1./hc->Integral());

        hCDF[s] = hc->GetCumulative();
        hCDF[s]->SetName(Form("hCDF_%s_%s", canvasName.Data(), histName[s].Data()));
        hCDF[s]->SetDirectory(0);
        hCDF[s]->SetStats(0);
        hCDF[s]->SetLineColor(sliceColor[s]);
        hCDF[s]->SetMarkerColor(sliceColor[s]);
        hCDF[s]->SetMarkerStyle(sliceMarker[s]);
    }

    hCDF[0]->SetTitle("");
    hCDF[0]->GetYaxis()->SetTitle("CDF");
    hCDF[0]->GetXaxis()->SetTitle(xTitle);
    hCDF[0]->GetYaxis()->SetRangeUser(0, 1);
    hCDF[0]->GetYaxis()->SetNdivisions(505);
    hCDF[0]->GetYaxis()->SetTitleFont(43);
    hCDF[0]->GetXaxis()->SetTitleFont(43);
    hCDF[0]->GetYaxis()->SetLabelFont(43);
    hCDF[0]->GetXaxis()->SetLabelFont(43);
    hCDF[0]->GetYaxis()->SetTitleSize(24);
    hCDF[0]->GetXaxis()->SetTitleSize(24);
    hCDF[0]->GetYaxis()->SetLabelSize(20);
    hCDF[0]->GetXaxis()->SetLabelSize(20);
    hCDF[0]->GetYaxis()->SetTitleOffset(1.3);
    hCDF[0]->GetXaxis()->SetTitleOffset(1.0);
    hCDF[0]->GetXaxis()->SetRangeUser(xDrawMin, xDrawMax);

    TLine *line = new TLine(xDrawMin, 0.5, xDrawMax, 0.5);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);

    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.05);

    TCanvas *c = new TCanvas(canvasName, canvasName, 800, 600);

    TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.3, 1.0, 1.0);
    pad1->SetLeftMargin(0.16); pad1->SetBottomMargin(0.15);
    pad1->Draw();

    TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.315);
    pad2->SetLeftMargin(0.16); pad2->SetBottomMargin(0.33); pad2->SetTopMargin(0.03);
    pad2->Draw();

    pad1->cd();
    if (islogY) pad1->SetLogy();
    h[0]->Draw("E0");
    for (int s = 1; s < nSlices; ++s) h[s]->Draw("E0 same");

    TLegend *leg = new TLegend(0.55, 0.68, 0.89, 0.88);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    for (int s = 0; s < nSlices; ++s) leg->AddEntry(h[s], etaLabel[s], "lp");
    leg->Draw();

    latex1->Draw();

    pad2->cd();
    gPad->SetTickx(0);
    hCDF[0]->Draw("hist P");
    for (int s = 1; s < nSlices; ++s) hCDF[s]->Draw("hist P same");
    line->Draw("same");

    c->cd();
    c->Modified();
    c->Update();

    if (outName.Length() > 0) {
        c->SaveAs(outDir + "/" + outName + ".pdf");
        latex1->SetTitle("#it{Private work (CMS data)}");
        c->Modified();
        c->Update();
        c->SaveAs(outDir + "/" + outName + "_bis.pdf");
        latex1->SetTitle("#scale[1.3]{#bf{CMS}}#it{Work in progress}");
    }

    return c;
}

void SignalMassResolution() {

    const std::string version = "19p8";

    const char *nameMass    = "METanalysis_TestPUppiMETCut_Eta2p4_9fp10_SignalMass_nominal";
    const char *nameMassATLAS = "METanalysis_TestPUppiMETCut_Eta2p4_9fp10_SignalMass_ATLASbkg";

    TFile *ofile = new TFile("../outputDisplay/SignalMassResolution.root", "RECREATE");

    // ================= Signal mass points =================
    std::vector<int> masses = {1100, 1400, 1600, 1800, 2000, 2200, 2400, 2600};

    std::vector<double> x, ex, y, ey;                 // nominal
    std::vector<double> xA, exA, yA, eyA;             // ATLAS

    for (int mass : masses) {
        std::string fname = "../output/Gluino_V19/Gluino_Run3_MET_madgraph_"
                          + std::to_string(mass) + "_V" + version + "_weighted.root";
        TFile *f = new TFile(fname.c_str(), "READ");
        if (!f || f->IsZombie()) { std::cerr << "Bad file for M=" << mass << std::endl; continue; }

        TH1D *h  = (TH1D*)f->Get(nameMass);
        TH1D *hA = (TH1D*)f->Get(nameMassATLAS);
        if (!h)  { std::cerr << "Missing " << nameMass      << " in " << f->GetName() << std::endl; f->Close(); continue; }
        if (!hA) { std::cerr << "Missing " << nameMassATLAS << " in " << f->GetName() << std::endl; }

        double mean  = h->GetMean();
        double rms   = h->GetRMS();      // width (sigma) de la distribution
        double meanE = h->GetMeanError();

        double ratio = (mass > 0) ? mean / mass : 0.;   // 1 = reco parfaite

        x.push_back(mass);
        ex.push_back(0.);
        y.push_back(mean);
        ey.push_back(rms);               // barre d'erreur = largeur de la masse reconstruite

        std::cout << "M = " << mass << " GeV [nominal] : mean = " << mean
                  << " , width (RMS) = " << rms
                  << " , meanErr = " << meanE
                  << " , mean/M = " << ratio
                  << " (" << (ratio - 1.) * 100. << " %)" << std::endl;

        if (hA) {
            double meanA  = hA->GetMean();
            double rmsA   = hA->GetRMS();
            double ratioA = (mass > 0) ? meanA / mass : 0.;

            xA.push_back(mass);
            exA.push_back(0.);
            yA.push_back(meanA);
            eyA.push_back(rmsA);

            std::cout << "M = " << mass << " GeV [ATLAS]   : mean = " << meanA
                      << " , width (RMS) = " << rmsA
                      << " , mean/M = " << ratioA
                      << " (" << (ratioA - 1.) * 100. << " %)" << std::endl;
        }

        f->Close();
    }

    ofile->cd();

    TGraphErrors *g  = new TGraphErrors(x.size(),  &x[0],  &y[0],  &ex[0],  &ey[0]);
    TGraphErrors *gA = new TGraphErrors(xA.size(), &xA[0], &yA[0], &exA[0], &eyA[0]);

    g->GetXaxis()->SetLimits(700., 2800.);
    g->GetYaxis()->SetRangeUser(700., 2800.);
    g->SetName("SignalMassResolution");
    gA->SetName("SignalMassResolution_ATLAS");

    // --- Canvas ---
    TCanvas *c = new TCanvas("c_SignalMassResolution", "Reconstructed mass vs true mass", 800, 600);
    c->SetGrid();
    c->SetLeftMargin(0.16); c->SetBottomMargin(0.16);

    // --- style nominal ---
    g->SetStats(0);
    g->SetTitle("");
    g->SetLineWidth(2);
    g->SetLineColor(kBlack);
    g->SetMarkerColor(kBlack);
    g->SetMarkerStyle(21);
    g->SetMarkerSize(1.2);
    g->SetFillColorAlpha(kRed+1, 0.50);   // bande de largeur

    g->GetXaxis()->SetTitle("Generated mass [GeV]");
    g->GetYaxis()->SetTitle("Reconstructed mass [GeV]");
    g->GetYaxis()->SetTitleSize(0.05);
    g->GetXaxis()->SetTitleSize(0.05);
    g->GetXaxis()->SetTitleOffset(1.1);
    g->GetYaxis()->SetTitleOffset(1.2);
    g->GetXaxis()->SetLabelSize(0.045);
    g->GetYaxis()->SetLabelSize(0.045);

    // --- style ATLAS ---
    gA->SetLineWidth(2);
    gA->SetLineColor(kBlue+1);
    gA->SetMarkerColor(kBlue+1);
    gA->SetMarkerStyle(20);
    gA->SetMarkerSize(1.2);
    gA->SetFillColorAlpha(kAzure+1, 0.40);   // bande de largeur

    // bande d'incertitude + points
    g->Draw("A3");        // bande remplie (largeur) nominal
    g->Draw("PX SAME");   // points + marqueurs nominal

    gA->Draw("3 SAME");   // bande remplie (largeur) ATLAS
    gA->Draw("PX SAME");  // points + marqueurs ATLAS

    // diagonale y = x pour reference
    TLine *diag = new TLine(700, 700, 2800, 2800);
    diag->SetLineStyle(2);
    diag->SetLineColor(kBlack);
    diag->Draw("SAME");

    TLegend *leg = new TLegend(0.2, 0.68, 0.55, 0.89);
    leg->SetBorderSize(0);
    leg->AddEntry(g,  "f^{KC} param. (mean #pm width)", "lpf");
    leg->AddEntry(gA, "f^{5p} param. (mean #pm width)", "lpf");
    //leg->AddEntry(diag, "y = x", "l");
    leg->Draw();

    TLatex *latex1 = new TLatex(0.16, 0.91, "#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.04);
    latex1->Draw();

    TLatex *latex2 = new TLatex(0.68, 0.91, "#scale[1.3]{#bf{m_{#tilde{g}}=2000 GeV}}");
    latex2->SetNDC();
    latex2->SetTextFont(42);
    latex2->SetTextSize(0.04);
    //latex2->Draw();

    c->SaveAs("../outputDisplay/SignalMassResolution.pdf");

    ofile->cd();
    g->Write();
    gA->Write();
    c->Write();

    latex1->SetTitle("#it{Private work (CMS simulation)}");
    c->Modified();
    c->Update();
    c->SaveAs("../outputDisplay/SignalMassResolution_bis.pdf");

    ofile->Close();

    return;
}


void etaSlices_bkg() {

    const TString inputPath =
        "/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/massSpectrum_bckgPrediction/DebugFit/Fits_HistForBkg_MC_9fp10_OldFit_Eta2p4_90_2.978400.root";

    const TString opHistName[nSlices] = { "forfit_p_eta10",
                                          "forfit_p_eta14",
                                          "forfit_p_eta20" };
    const TString ihHistName[nSlices] = { "proj_ih_eta10",
                                          "proj_ih_eta14",
                                          "proj_ih_eta20" };

    const double opDrawMax = 150.;  // x-axis display range for 1/p
    const double ihDrawMax = 5.5;   // x-axis display range for Ih

    TFile *ifile = TFile::Open(inputPath, "READ");
    if (!ifile || ifile->IsZombie()) {
        Error("etaSlices_bkg", "cannot open %s", inputPath.Data());
        return;
    }

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    drawEtaSlices(ifile, opHistName, "10^{4}/p [GeV^{-1}]",
                  0., opDrawMax, 50., "c_oP", "oP_bkg_etaSlices");

    drawEtaSlices(ifile, ihHistName, "I_{h} [MeV/cm]",
                  2.97, ihDrawMax, 50., "c_Ih", "Ih_bkg_etaSlices");

    ifile->Close();   // safe: the histograms are cloned and detached

    return;
}


namespace SelEff {

// ---------------------------------------------------------------------------
// Recherche recursive d'un TH1 dans le fichier (les histos HSCP sont souvent
// ranges dans un TDirectory de type "HSCParticleAnalyzer/BaseName/...").
// suffixMatch : accepte un nom prefixe (ex. "METanalysis_Nosel_P" pour "Nosel_P").
// ---------------------------------------------------------------------------
TH1 *FindHisto(TDirectory *dir, const TString &name, bool suffixMatch) {
   if (!dir) return nullptr;

   if (!suffixMatch) {
      TObject *obj = dir->Get(name); // gere aussi les chemins "sousdir/nom"
      if (obj && obj->InheritsFrom(TH1::Class())) return static_cast<TH1 *>(obj);
   }

   TIter next(dir->GetListOfKeys());
   while (TKey *key = static_cast<TKey *>(next())) {
      TClass *cl = TClass::GetClass(key->GetClassName());
      if (!cl) continue;
      if (cl->InheritsFrom(TDirectory::Class())) {
         TDirectory *sub = static_cast<TDirectory *>(key->ReadObj());
         if (TH1 *h = FindHisto(sub, name, suffixMatch)) return h;
      } else if (suffixMatch && cl->InheritsFrom(TH1::Class())) {
         TString kn(key->GetName());
         if (kn == name || kn.EndsWith(name)) return static_cast<TH1 *>(key->ReadObj());
      }
   }
   return nullptr;
}

TH1 *GetHisto(TFile *f, const TString &name) {
   TH1 *h = FindHisto(f, name, false);
   if (!h) h = FindHisto(f, name, true);
   return h;
}

// ---------------------------------------------------------------------------
// Efficacite + erreur binomiale calculee sur les entrees effectives
// (robuste si les histogrammes sont ponderes : neff = (sum w)^2 / sum w^2).
// ---------------------------------------------------------------------------
struct Eff {
   double val = 0.;
   double err = 0.;
   bool ok = false;
};

Eff ComputeEff(const TH1 *hSel, const TH1 *hTot, bool withOverflow) {
   Eff e;
   if (!hSel || !hTot) return e;

   const int lo = withOverflow ? 0 : 1;
   const int hiSel = withOverflow ? hSel->GetNbinsX() + 1 : hSel->GetNbinsX();
   const int hiTot = withOverflow ? hTot->GetNbinsX() + 1 : hTot->GetNbinsX();

   const double nSel = hSel->Integral(lo, hiSel);
   const double nTot = hTot->Integral(lo, hiTot);
   if (nTot <= 0.) return e;

   double sumw2 = 0.;
   for (int i = lo; i <= hiTot; ++i) {
      const double be = hTot->GetBinError(i);
      sumw2 += be * be;
   }
   const double nEff = (sumw2 > 0.) ? nTot * nTot / sumw2 : nTot;

   e.val = nSel / nTot;
   const double var = e.val * (1. - e.val) / nEff;
   e.err = (var > 0.) ? std::sqrt(var) : 0.;
   e.ok = true;
   return e;
}

struct EtaCat {
   const char *tag;
   const char *label;
};

struct Variant {
   const char *tag;
   const char *label;
   int color;
   int marker;
};

} // namespace SelEff

// ---------------------------------------------------------------------------
void CompareSelEff(TString inputDir = ".", TString outDir = "SelEff",
                   TString noselName = "Nosel_P", bool withOverflow = true,
                   bool logY = false) {
   using namespace SelEff;

   const std::vector<int> masses = {1100, 1200, 1300, 1400, 1600,
                                    1800, 2000, 2200, 2400, 2600};
   const char *fileFmt = "Gluino_Run3_MET_madgraph_%d_V19p12.root";
   const TString base = "METanalysis_TestPUppiMETCut_";
   const TString hSuffix = "_PUppiMET_VS_PseudoMET";

   const std::vector<EtaCat> etaCats = {{"Eta2p4", "|#eta| < 2.4"},
                                        {"Eta1", "|#eta| < 1.0"},
                                        {"Eta1_2p4", "1.0 < |#eta| < 2.4"}};

   const std::vector<Variant> variants = {
       {"", "Nominal", kBlack, 20},
       {"SigmaPtoverPt_0p5_", "#sigma_{p_{T}}/p_{T} < 0.5", kAzure + 2, 21},
       {"EoP_0p1_", "E/p < 0.1", kRed + 1, 22},
       {"SigmaPtoverPt_0p5_EoP_0p1_", "#sigma_{p_{T}}/p_{T} < 0.5 & E/p < 0.1",
        kGreen + 2, 23}};

   const size_t nEta = etaCats.size();
   const size_t nVar = variants.size();

   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   gSystem->mkdir(outDir, kTRUE);

   // --- creation des graphes -----------------------------------------------
   std::vector<std::vector<TGraphErrors *>> gr(
       nEta, std::vector<TGraphErrors *>(nVar, nullptr));
   for (size_t ie = 0; ie < nEta; ++ie) {
      for (size_t iv = 0; iv < nVar; ++iv) {
         TGraphErrors *g = new TGraphErrors();
         g->SetName(Form("eff_%s_%s", etaCats[ie].tag,
                         (strlen(variants[iv].tag) ? variants[iv].tag : "nominal")));
         g->SetTitle(variants[iv].label);
         g->SetLineColor(variants[iv].color);
         g->SetMarkerColor(variants[iv].color);
         g->SetMarkerStyle(variants[iv].marker);
         g->SetMarkerSize(1.2);
         g->SetLineWidth(2);
         gr[ie][iv] = g;
      }
   }

   // --- boucle sur les points de masse -------------------------------------
   for (int m : masses) {
      const TString fpath =
          TString::Format("%s/%s", inputDir.Data(), TString::Format(fileFmt, m).Data());
      TFile *f = TFile::Open(fpath, "READ");
      if (!f || f->IsZombie()) {
         std::cerr << "[CompareSelEff] fichier manquant : " << fpath << std::endl;
         if (f) delete f;
         continue;
      }

      TH1 *hTot = GetHisto(f, noselName);
      if (!hTot) {
         std::cerr << "[CompareSelEff] denominateur '" << noselName
                   << "' introuvable dans " << fpath << std::endl;
         f->Close();
         delete f;
         continue;
      }

      std::cout << "\n=== m_gluino = " << m << " GeV  (N_Nosel = "
                << hTot->Integral(withOverflow ? 0 : 1,
                                  withOverflow ? hTot->GetNbinsX() + 1 : hTot->GetNbinsX())
                << ") ===" << std::endl;

      for (size_t ie = 0; ie < nEta; ++ie) {
         for (size_t iv = 0; iv < nVar; ++iv) {
            const TString selName =
                base + variants[iv].tag + etaCats[ie].tag + hSuffix;
            TH1 *hSel = GetHisto(f, selName);
            if (!hSel) {
               std::cerr << "  [manquant] " << selName << std::endl;
               continue;
            }
            const Eff e = ComputeEff(hSel, hTot, withOverflow);
            if (!e.ok) continue;

            const int n = gr[ie][iv]->GetN();

            gr[ie][iv]->SetPoint(n, m, e.val*0.855);
            gr[ie][iv]->SetPointError(n, 0., e.err*0.855);

            printf("  %-12s %-30s eff = %7.4f +/- %6.4f\n", etaCats[ie].tag,
                   (strlen(variants[iv].tag) ? variants[iv].tag : "nominal"), e.val,
                   e.err);
         }
      }

      f->Close();
      delete f;
   }

   // --- trace : un canvas par categorie en eta ------------------------------
   TFile *fout = TFile::Open(outDir + "/SelEff_graphs.root", "RECREATE");

   for (size_t ie = 0; ie < nEta; ++ie) {
      TCanvas *c = new TCanvas(Form("c_%s", etaCats[ie].tag),
                               Form("Efficacite %s", etaCats[ie].tag), 900, 600);
      c->SetTicks(1, 1);
      c->SetGridy();
      c->SetLeftMargin(0.16);
      c->SetBottomMargin(0.16);
      if (logY) c->SetLogy();

      TMultiGraph *mg = new TMultiGraph();
      double ymax = 0.;
      for (size_t iv = 0; iv < nVar; ++iv) {
         if (gr[ie][iv]->GetN() == 0) continue;
         mg->Add(gr[ie][iv], "PL");
         for (int i = 0; i < gr[ie][iv]->GetN(); ++i)
            ymax = std::max(ymax, gr[ie][iv]->GetY()[i] + gr[ie][iv]->GetEY()[i]);
      }
      if (mg->GetListOfGraphs() == nullptr) {
         std::cerr << "[CompareSelEff] aucun point pour " << etaCats[ie].tag << std::endl;
         continue;
      }

      mg->Draw("AP");
      mg->GetXaxis()->SetTitle("m_{#tilde{g}} [GeV]");
      mg->GetYaxis()->SetTitle("Signal acceptance");
        mg->GetYaxis()->SetTitleSize(0.06);
        mg->GetXaxis()->SetTitleSize(0.06);
        mg->GetXaxis()->SetTitleOffset(0.9);
        mg->GetYaxis()->SetTitleOffset(1);
        mg->GetXaxis()->SetLabelSize(0.05);
        mg->GetYaxis()->SetLabelSize(0.05);
        TLatex *latex1 = new TLatex(0.16, 0.91, "#it{Private work (CMS simulation)}");
        latex1->SetNDC();
        latex1->SetTextFont(42);
        latex1->SetTextSize(0.04);
        latex1->Draw("same");
      mg->GetXaxis()->SetLimits(masses.front() - 100, masses.back() + 100);
      mg->SetMinimum(ymax*0.7);
      mg->SetMaximum(ymax * 1.3);
      c->Modified();

      TLegend *leg = new TLegend(0.18, 0.68, 0.62, 0.88);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.030);
      for (size_t iv = 0; iv < nVar; ++iv)
         if (gr[ie][iv]->GetN() > 0) leg->AddEntry(gr[ie][iv], variants[iv].label, "lp");
      leg->Draw();

      TLatex tex;
      tex.SetNDC();
      tex.SetTextFont(42);
      tex.SetTextSize(0.038);
      tex.DrawLatex(0.75, 0.91, Form("#bf{%s}",etaCats[ie].label));

      c->Update();
      c->SaveAs(Form("%s/SelEff_%s.pdf", outDir.Data(), etaCats[ie].tag));

      if (fout && !fout->IsZombie()) {
         fout->cd();
         for (size_t iv = 0; iv < nVar; ++iv)
            if (gr[ie][iv]->GetN() > 0) gr[ie][iv]->Write();
         c->Write();
      }
   }

   if (fout) {
      fout->Close();
      delete fout;
   }
   std::cout << "\n[CompareSelEff] sorties ecrites dans " << outDir << std::endl;
}

// ---------------------------------------------------------------------------
// A ajouter dans CombineHistos.C, apres CompareSelEff().
// Un PDF multi-pages par categorie en eta ; une page par point de masse ;
// 4 histogrammes (les 4 variantes de selection) superposes par page.
// Histogramme utilise : <selLabel>_9fp10_SignalMass_nominal
// ---------------------------------------------------------------------------
void PlotMassSpectraPerEta(TString inputDir = ".", TString outDir = "SelEff",
                           bool normalize = false, bool logY = true,
                           int rebin = 1) {
   using namespace SelEff;

   const std::vector<int> masses = {1100, 1200, 1300, 1400, 1600,
                                    1800, 2000, 2200, 2400, 2600};
   const char *fileFmt = "Gluino_Run3_MET_madgraph_%d_V19p12.root";
   const TString base = "METanalysis_TestPUppiMETCut_";
   const TString hSuffix = "_9fp10_SignalMass_nominal";

   const std::vector<EtaCat> etaCats = {{"Eta2p4", "|#eta| < 2.4"},
                                        {"Eta1", "|#eta| < 1.0"},
                                        {"Eta1_2p4", "1.0 < |#eta| < 2.4"}};

   const std::vector<Variant> variants = {
       {"", "Nominal", kBlack, 20},
       {"SigmaPtoverPt_0p5_", "#sigma_{p_{T}}/p_{T} < 0.5", kAzure + 2, 21},
       {"EoP_0p1_", "E/p < 0.1", kRed + 1, 22},
       {"SigmaPtoverPt_0p5_EoP_0p1_", "#sigma_{p_{T}}/p_{T} < 0.5 & E/p < 0.1",
        kGreen + 2, 23}};

   const size_t nEta = etaCats.size();
   const size_t nVar = variants.size();

   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   gSystem->mkdir(outDir, kTRUE);

   for (size_t ie = 0; ie < nEta; ++ie) {
      const TString pdfName =
          Form("%s/MassSpectra_%s.pdf", outDir.Data(), etaCats[ie].tag);

      TCanvas *c = new TCanvas(Form("cMass_%s", etaCats[ie].tag), "", 800, 600);
      c->SetTicks(1, 1);
      c->SetLeftMargin(0.12);
      c->SetBottomMargin(0.12);
      if (logY) c->SetLogy();
      c->Print(pdfName + "["); // ouverture du PDF multi-pages

      for (int m : masses) {
         const TString fpath = TString::Format(
             "%s/%s", inputDir.Data(), TString::Format(fileFmt, m).Data());
         TFile *f = TFile::Open(fpath, "READ");
         if (!f || f->IsZombie()) {
            std::cerr << "[MassSpectra] fichier manquant : " << fpath << std::endl;
            if (f) delete f;
            continue;
         }

         // --- recuperation + clonage (les histos doivent survivre au Close) ---
         std::vector<TH1 *> hs(nVar, nullptr);
         for (size_t iv = 0; iv < nVar; ++iv) {
            const TString hName =
                base + variants[iv].tag + etaCats[ie].tag + hSuffix;
            TH1 *h0 = GetHisto(f, hName);
            if (!h0) {
               std::cerr << "  [manquant] " << hName << " (m = " << m << ")"
                         << std::endl;
               continue;
            }
            TH1 *h = static_cast<TH1 *>(h0->Clone(
                Form("mass_%d_%s_%zu", m, etaCats[ie].tag, iv)));
            h->SetDirectory(nullptr);
            h = rebinHisto(h);
            h->SetDirectory(nullptr);
            if (rebin > 1) h->Rebin(rebin);
            const double integ = h->Integral(0, h->GetNbinsX() + 1);
            if (normalize && integ > 0.) h->Scale(1. / integ);
            h->SetLineColor(variants[iv].color);
            h->SetMarkerColor(variants[iv].color);
            h->SetMarkerStyle(variants[iv].marker);
            h->SetMarkerSize(0.9);
            h->SetLineWidth(2);
            h->SetStats(0);
            hs[iv] = h;
         }
         f->Close();
         delete f;

         double ymax = 0.;
         for (size_t iv = 0; iv < nVar; ++iv)
            if (hs[iv]) ymax = std::max(ymax, hs[iv]->GetMaximum());
         if (ymax <= 0.) {
            for (auto *h : hs) delete h;
            continue;
         }

         // --- trace de la page ------------------------------------------------
         c->cd();
         c->Clear();
         c->SetTicks(1, 1);
         c->SetLeftMargin(0.12);
         c->SetBottomMargin(0.12);
         if (logY) c->SetLogy();

         TLegend *leg = new TLegend(0.45, 0.13, 0.8, 0.4);
         leg->SetBorderSize(0);
         leg->SetFillStyle(0);
         leg->SetTextSize(0.030);

         bool first = true;
         for (size_t iv = 0; iv < nVar; ++iv) {
            if (!hs[iv]) continue;
            if (first) {
               hs[iv]->GetXaxis()->SetTitle("m [GeV]");
               hs[iv]->GetXaxis()->SetRangeUser(0,4000);
               hs[iv]->GetYaxis()->SetTitle(normalize ? "a.u." : "Events");
               hs[iv]->GetXaxis()->SetTitleSize(0.06);
               hs[iv]->GetYaxis()->SetTitleSize(0.06);
               hs[iv]->GetXaxis()->SetTitleOffset(0.9);
               hs[iv]->GetYaxis()->SetTitleOffset(0.9);
               hs[iv]->GetXaxis()->SetLabelSize(0.05);
               hs[iv]->GetYaxis()->SetLabelSize(0.05);
               hs[iv]->SetMinimum(logY ? (normalize ? 1e-5 : 0.5) : 0.);
               hs[iv]->SetMaximum(logY ? ymax * 20. : ymax * 1.45);

               hs[iv]->Draw("HIST E");
               first = false;
            } else {
               hs[iv]->Draw("HIST E SAME");
            }
            leg->AddEntry(hs[iv],
                          Form("%s  (N = %.0f)", variants[iv].label,
                               hs[iv]->Integral(0, hs[iv]->GetNbinsX() + 1)),
                          "l");
         }
         leg->Draw();

         TLatex tex;
         tex.SetNDC();
         tex.SetTextFont(42);
         tex.SetTextSize(0.04);
         tex.DrawLatex(0.12, 0.91, "#it{Private work (CMS simulation)}");
         tex.SetTextSize(0.038);
         tex.DrawLatex(0.62, 0.91,
                       Form("#bf{m_{#tilde{g}} = %d GeV, %s}", m,
                            etaCats[ie].label));

         c->Modified();
         c->Update();
         c->Print(pdfName); // une page

         delete leg;
         for (auto *h : hs) delete h;
      }

      c->Print(pdfName + "]"); // fermeture du PDF
      delete c;
      std::cout << "[MassSpectra] " << pdfName << " ecrit." << std::endl;
   }
}

struct LinFitResult {
    double p0 = 0., p0Err = 0.;   // ordonnee a l'origine
    double p1 = 0., p1Err = 0.;   // pente
    double cov01 = 0.;            // covariance p0-p1 (utile pour propager l'erreur)
    double chi2 = 0.;
    int    ndf  = 0;
    double prob = 0.;
    int    nPoints = 0;
    bool   ok = false;
};

LinFitResult RatioLinearFit(const TH1* hNum_,
                            const TH1* hDen_,
                            const char* outName,
                            bool   normalise = true,
                            double xmin = 0., double xmax = -1.,
                            double minDenContent = 0.,
                            TH1D** hRatioOut = nullptr) {
    LinFitResult res;

    if (!hNum_ || !hDen_) {
        std::cerr << "RatioLinearFit(" << outName << "): histogramme nul" << std::endl;
        return res;
    }

    // --- Compatibilite de binning ---
    const int nb = hNum_->GetNbinsX();
    if (nb != hDen_->GetNbinsX() ||
        std::fabs(hNum_->GetXaxis()->GetXmin() - hDen_->GetXaxis()->GetXmin()) > 1e-9 ||
        std::fabs(hNum_->GetXaxis()->GetXmax() - hDen_->GetXaxis()->GetXmax()) > 1e-9) {
        std::cerr << "RatioLinearFit(" << outName << "): binnings incompatibles ("
                  << nb << " vs " << hDen_->GetNbinsX() << " bins)" << std::endl;
        return res;
    }

    // --- Copies de travail ---
    std::unique_ptr<TH1D> hNum(static_cast<TH1D*>(hNum_->Clone(Form("%s_num", outName))));
    std::unique_ptr<TH1D> hDen(static_cast<TH1D*>(hDen_->Clone(Form("%s_den", outName))));
    hNum->SetDirectory(nullptr);
    hDen->SetDirectory(nullptr);

    if (normalise) {
        const double iN = hNum->Integral(0, nb + 1);
        const double iD = hDen->Integral(0, nb + 1);
        if (iN <= 0. || iD <= 0.) {
            std::cerr << "RatioLinearFit(" << outName << "): integrale nulle ou negative"
                      << std::endl;
            return res;
        }
        hNum->Scale(1. / iN);
        hDen->Scale(1. / iD);
    }

    // --- Rapport bin a bin, en ne gardant que les bins exploitables ---
    std::unique_ptr<TH1D> hRatio(static_cast<TH1D*>(hNum->Clone(outName)));
    hRatio->SetDirectory(nullptr);
    hRatio->Reset("ICESM");
    hRatio->SetTitle(Form("%s;%s;ratio", outName, hNum_->GetXaxis()->GetTitle()));

    std::vector<double> vx, vy, vex, vey;
    vx.reserve(nb); vy.reserve(nb); vex.reserve(nb); vey.reserve(nb);

    for (int i = 1; i <= nb; ++i) {
        const double n  = hNum->GetBinContent(i);
        const double d  = hDen->GetBinContent(i);
        const double en = hNum->GetBinError(i);
        const double ed = hDen->GetBinError(i);

        if (d <= 0. || d < minDenContent) continue;   // dividende impossible
        if (n <= 0.)                      continue;   // bin vide cote numerateur

        const double r = n / d;
        // propagation lineaire, numerateur et denominateur supposes independants
        const double er = r * std::sqrt(std::pow(en / n, 2) + std::pow(ed / d, 2));

        hRatio->SetBinContent(i, r);
        hRatio->SetBinError(i, er);

        const double xc = hRatio->GetBinCenter(i);
        if (xmax > xmin && (xc < xmin || xc > xmax)) continue;  // hors plage de fit

        vx.push_back(xc);
        vy.push_back(r);
        vex.push_back(0.);
        vey.push_back(er > 0. ? er : 1e-9);
    }

    res.nPoints = static_cast<int>(vx.size());
    if (res.nPoints < 3) {
        std::cerr << "RatioLinearFit(" << outName << "): seulement " << res.nPoints
                  << " points exploitables, fit abandonne" << std::endl;
        if (hRatioOut) *hRatioOut = static_cast<TH1D*>(hRatio.release());
        return res;
    }

    // --- Fit lineaire ---
    TGraphErrors g(res.nPoints, vx.data(), vy.data(), vex.data(), vey.data());
    g.SetName(Form("%s_gr", outName));

    const double fitLo = (xmax > xmin) ? xmin : vx.front();
    const double fitHi = (xmax > xmin) ? xmax : vx.back();

    TF1 fLin(Form("%s_fLin", outName), "[0]+[1]*x", fitLo, fitHi);
    fLin.SetParameters(1., 0.);
    fLin.SetParNames("p0", "p1");

    TFitResultPtr fr = g.Fit(&fLin, "QRS");
    if (fr.Get() == nullptr || !fr->IsValid()) {
        std::cerr << "RatioLinearFit(" << outName << "): le fit n'a pas converge" << std::endl;
        if (hRatioOut) *hRatioOut = static_cast<TH1D*>(hRatio.release());
        return res;
    }

    res.p0    = fr->Parameter(0);  res.p0Err = fr->ParError(0);
    res.p1    = fr->Parameter(1);  res.p1Err = fr->ParError(1);
    res.cov01 = fr->CovMatrix(0, 1);
    res.chi2  = fr->Chi2();
    res.ndf   = fr->Ndf();
    res.prob  = (res.ndf > 0) ? TMath::Prob(res.chi2, res.ndf) : -1.;
    res.ok    = true;

    // la fonction reste attachee au rapport pour le dessin
    hRatio->GetListOfFunctions()->Add(fLin.Clone());

    std::cout << "  [" << outName << "]  p0 = " << res.p0 << " +/- " << res.p0Err
              << " | p1 = " << res.p1 << " +/- " << res.p1Err
              << " | chi2/ndf = " << res.chi2 << "/" << res.ndf
              << " (p = " << res.prob << ")"
              << " | " << res.nPoints << " pts";
    if (res.p1Err > 0.)
        std::cout << " | pente a " << std::fabs(res.p1) / res.p1Err << " sigma de 0";
    std::cout << std::endl;

    if (hRatioOut) *hRatioOut = static_cast<TH1D*>(hRatio.release());
    return res;
}

// Canvas a deux pads : histos superposes en haut, rapport + fit en bas.
// Le TF1 est deja attache a hRatio par RatioLinearFit, il se dessine tout seul.
TCanvas* MakeRatioCanvas(const TH1* hNum_, const TH1* hDen_, const TH1D* hRatio_,
                         const LinFitResult& fit,
                         const char* cname,
                         const char* legNum, const char* legDen,
                         const std::string& XaxisTitle = "",
                         float Xmin = 0., float Xmax = -1.,
                         float RatioMin = 0.5, float RatioMax = 1.5,
                         bool normalise = true) {
    if (!hNum_ || !hDen_ || !hRatio_) return nullptr;

    gStyle->SetOptStat(0);

    TCanvas* c_new = new TCanvas(cname, cname, 800, 600);

    TPad* pad1 = new TPad(Form("%s_pad1", cname), "pad1", 0.0, 0.30, 1.0, 1.0);
    pad1->SetLeftMargin(0.16);
    pad1->SetBottomMargin(0.15);
    pad1->Draw();

    TPad* pad2 = new TPad(Form("%s_pad2", cname), "pad2", 0.0, 0.0, 1.0, 0.315);
    pad2->SetLeftMargin(0.16);
    pad2->SetBottomMargin(0.33);
    pad2->SetTopMargin(0.02);
    pad2->Draw();

    // bornes en x : auto si non fournies
    if (Xmax <= Xmin) {
        Xmin = hRatio_->GetXaxis()->GetXmin();
        Xmax = hRatio_->GetXaxis()->GetXmax();
    }

    // ---------------- Pad du haut ----------------
    pad1->cd();

    TH1D* hN = (TH1D*)hNum_->Clone(Form("%s_drawN", cname));
    TH1D* hD = (TH1D*)hDen_->Clone(Form("%s_drawD", cname));
    hN->SetDirectory(nullptr);  hN->SetBit(kCanDelete);
    hD->SetDirectory(nullptr);  hD->SetBit(kCanDelete);

    if (normalise) {   // meme normalisation que dans RatioLinearFit
        const double iN = hN->Integral(0, hN->GetNbinsX() + 1);
        const double iD = hD->Integral(0, hD->GetNbinsX() + 1);
        if (iN > 0.) hN->Scale(1. / iN);
        if (iD > 0.) hD->Scale(1. / iD);
    }

    hN->SetTitle("");
    hN->SetLineColor(kBlue);  hN->SetMarkerColor(kBlue);  hN->SetMarkerStyle(8);
    hD->SetLineColor(kRed);    hD->SetMarkerColor(kRed);    hD->SetMarkerStyle(24);

    hN->GetXaxis()->SetTitleFont(43);
    hN->GetXaxis()->SetLabelFont(43);
    hN->GetXaxis()->SetTitleSize(24);
    hN->GetXaxis()->SetLabelSize(20);
    hN->GetYaxis()->SetTitleFont(43);
    hN->GetYaxis()->SetLabelFont(43);
    hN->GetYaxis()->SetTitleSize(24);     // px
    hN->GetYaxis()->SetLabelSize(20);
    hN->GetYaxis()->SetTitleOffset(1.8);
    hN->GetXaxis()->SetTitleOffset(0.9);
    hN->GetXaxis()->SetTitle(XaxisTitle.c_str());
    hN->GetYaxis()->SetTitle(normalise ? "Events (normalised)" : "Events");
    hN->GetXaxis()->SetRangeUser(Xmin, Xmax);
    hN->SetMaximum(5 * std::max(hN->GetMaximum(), hD->GetMaximum()));

    hN->Draw("E0");
    hD->Draw("E0 same");

    TLegend* leg = new TLegend(0.75, 0.7, 0.9, 0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(43);
    leg->SetTextSize(20);
    leg->AddEntry(hN, legNum, "lep");
    leg->AddEntry(hD, legDen, "lep");
    leg->SetBit(kCanDelete);
    leg->Draw();

    pad1->SetLogy();

    // ---------------- Pad du bas ----------------
    pad2->cd();

    TH1D* hR = (TH1D*)hRatio_->Clone(Form("%s_drawR", cname));
    hR->SetDirectory(nullptr);
    hR->GetListOfFunctions()->Delete();   // on jette le TF1 clone par Clone()
    hR->SetBit(kCanDelete);

    hR->SetTitle("");
    hR->SetLineColor(kBlack);
    hR->SetMarkerColor(kBlack);
    hR->SetMarkerStyle(8);

    hR->GetYaxis()->SetTitle("in CR/ in VR");
    hR->GetXaxis()->SetTitle(XaxisTitle.c_str());
    hR->GetYaxis()->SetNdivisions(505);
    hR->GetYaxis()->SetTitleFont(43);   hR->GetXaxis()->SetTitleFont(43);
    hR->GetYaxis()->SetLabelFont(43);   hR->GetXaxis()->SetLabelFont(43);
    hR->GetYaxis()->SetTitleSize(24);   hR->GetXaxis()->SetTitleSize(24);
    hR->GetYaxis()->SetLabelSize(20);   hR->GetXaxis()->SetLabelSize(20);
    hR->GetYaxis()->SetTitleOffset(1.6);
    hR->GetXaxis()->SetTitleOffset(1);
    hR->GetYaxis()->SetRangeUser(RatioMin, RatioMax);

    gPad->SetTickx(0);
    hR->Draw("E0");
    hR->GetXaxis()->SetRangeUser(Xmin, Xmax);

    TF1* fDraw = new TF1(Form("%s_fDraw", cname), "[0]+[1]*x", Xmin, Xmax);
    fDraw->SetParameters(fit.p0, fit.p1);
    fDraw->SetLineColor(kRed);
    fDraw->SetLineWidth(2);
    fDraw->SetBit(kCanDelete);
    fDraw->Draw("same");

    TLine* line = new TLine(hRatio_->GetXaxis()->GetXmin(), 1., Xmax, 1.);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    line->SetBit(kCanDelete);
    line->Draw("same");

    c_new->cd();
    c_new->Update();

    pad1->cd();
    TLatex* latex1 = new TLatex(0.16, 0.91, "#it{Private work (CMS data)}");
    latex1->SetNDC();
    latex1->SetTextFont(43);
    latex1->SetTextSize(24);
    latex1->SetBit(kCanDelete);
    latex1->Draw();
    pad1->Modified();

    std::cout << "Canvas " << cname << " dessine : " << legNum << " / " << legDen << std::endl;
    return c_new;
}

struct Pair { const char* num; const char* den; const char* out; const char *name;};

std::map<std::string, LinFitResult> FitVRRatios(const std::string& filename,
            const std::string& suffix = "",
            bool normalise = true,
            bool writeOut  = true) {
    std::map<std::string, LinFitResult> results;

    if (gROOT->GetListOfFiles()->FindObject(filename.c_str())) {
        std::cerr << "FitVRRatios: " << filename << " est deja ouvert" << std::endl;
        return results;
    }

    std::unique_ptr<TFile> f(TFile::Open(filename.c_str(), writeOut ? "UPDATE" : "READ"));
    if (!f || f->IsZombie()) {
        std::cerr << "FitVRRatios: impossible d'ouvrir " << filename << std::endl;
        return results;
    }

    const std::vector<Pair> pairs = {
        { "ih_eta_mean", "ih_VR",  "ratio_ih_meanOverVR", "I_{h} [MeV/cm]" },
        { "oP_eta_mean", "eta_VR", "ratio_oP_meanOverVR", "10^{4}/p [GeV^{-1}]" }
    };

    for (const auto& p : pairs) {
        TH1* hN = dynamic_cast<TH1*>(f->Get((std::string(p.num)).c_str()));
        TH1* hD = dynamic_cast<TH1*>(f->Get((std::string(p.den)).c_str()));
        if (!hN || !hD) {
            std::cerr << "FitVRRatios: introuvable ->"
                      << (hN ? "" : (std::string(" ") + p.num))
                      << (hD ? "" : (std::string(" ") + p.den)) << std::endl;
            continue;
        }

        const std::string outName = std::string(p.out) + suffix;

        TH1D* hRatio = nullptr;
        LinFitResult r = RatioLinearFit(hN, hD, outName.c_str(),
                                        normalise, 0., -1., 0., &hRatio);
        results[outName] = r;

        if (hRatio) {
            TCanvas* c = MakeRatioCanvas(hN, hD, hRatio, r,
                                         ("c_" + outName).c_str(),
                                         "in CR", "in VR",
                                         p.name, 0., (std::strcmp(p.den, "ih_VR") == 0) ? 6 : 160, 0.5, 1.5, normalise);
            if (writeOut && c) {
                f->cd();
                c->Write("", TObject::kOverwrite);
                c->SaveAs(("c_" + outName + ".pdf").c_str());   // optionnel
            }
        }
    }
    
    return results;
}

// ---------------------------------------------------------------------------
// Run a ROOT-generated .C macro (whose function name == the file stem) and
// return the LAST canvas created. This avoids hard-coding the internal canvas
// name, which differs between PseudoCaloMET / PUppiMET / ... files.
// ---------------------------------------------------------------------------
static TCanvas* RunAndGetLastCanvas(const char *macroPath) {
    int nBefore = gROOT->GetListOfCanvases()->GetSize();
    gROOT->ProcessLine(Form(".x %s", macroPath));
    TList *cList = (TList*)gROOT->GetListOfCanvases();
    if (cList->GetSize() <= nBefore) {
        std::cerr << "Error: no new canvas created by " << macroPath << std::endl;
        return nullptr;
    }
    return (TCanvas*)cList->Last();
}

// Find the first TH1 in a pad whose name contains 'needle'.
static TH1* FindHistInPad(TPad *pad, const char *needle) {
    if (!pad) return nullptr;
    TIter next(pad->GetListOfPrimitives());
    TObject *obj;
    while ((obj = next())) {
        if (obj->InheritsFrom(TH1::Class()) &&
            TString(obj->GetName()).Contains(needle))
            return (TH1*)obj;
    }
    return nullptr;
}

// Find the first TH1 in a canvas (single-pad signal file), scanning sub-pads too.
static TH1* FindHistInCanvas(TCanvas *c, const char *needle) {
    if (!c) return nullptr;
    TH1 *h = FindHistInPad((TPad*)c, needle);
    if (h) return h;
    TIter next(c->GetListOfPrimitives());
    TObject *obj;
    while ((obj = next())) {
        if (obj->InheritsFrom(TPad::Class())) {
            h = FindHistInPad((TPad*)obj, needle);
            if (h) return h;
        }
    }
    return nullptr;
}

// Find the top pad: the TPad named "pad1", else the first TPad found.
static TPad* FindTopPad(TCanvas *c) {
    if (!c) return nullptr;
    TPad *p = (TPad*)c->GetListOfPrimitives()->FindObject("pad1");
    if (p) return p;
    TIter next(c->GetListOfPrimitives());
    TObject *obj;
    while ((obj = next()))
        if (obj->InheritsFrom(TPad::Class())) return (TPad*)obj;
    return nullptr;
}

// Output stem from the ratio macro path (strip dir + ".C").
static TString OutStem(const char *ratioMacro) {
    TString s(ratioMacro);
    Ssiz_t slash = s.Last('/');
    if (slash != kNPOS) s.Remove(0, slash + 1);
    if (s.EndsWith(".C")) s.Remove(s.Length() - 2);
    return s;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
void OverlaySignalOnRatio(const char *signalMacro,
                          const char *ratioMacro,
                          const char *signalLabel = "Signal",
                          const char *outDir      = "TriggEff") {

    gErrorIgnoreLevel = kWarning;

    // 1) SIGNAL macro -> histogram whose name contains "SIGNAL"
    TCanvas *cSig = RunAndGetLastCanvas(signalMacro);
    if (!cSig) return;
    TH1 *hSig = FindHistInCanvas(cSig, "SIGNAL");
    if (!hSig) {
        std::cerr << "Error: no 'SIGNAL' histogram in " << signalMacro << std::endl;
        return;
    }
    TH1 *hSigClone = (TH1*)hSig->Clone("hSig_overlay");
    hSigClone->SetDirectory(0);

    // 2) RATIO macro -> canvas + top pad
    TCanvas *cRatio = RunAndGetLastCanvas(ratioMacro);
    if (!cRatio) return;
    TPad *pad1 = FindTopPad(cRatio);
    if (!pad1) {
        std::cerr << "Error: top pad not found in ratio canvas." << std::endl;
        return;
    }

    // sanity: warn if axes differ (PUppi vs PseudoCalo)
    TH1 *hData = FindHistInPad(pad1, "DATA");
    if (hData) {
        TString sx = hSigClone->GetXaxis()->GetTitle();
        TString dx = hData->GetXaxis()->GetTitle();
        if (sx != dx)
            std::cerr << "Note: signal x-axis '" << sx << "' != data x-axis '"
                      << dx << "'. Overlaying anyway." << std::endl;
    }

    // 3) style + draw signal on top pad
    hSigClone->SetLineColor(kGreen + 2);
    hSigClone->SetMarkerColor(kGreen + 2);
    hSigClone->SetMarkerStyle(21);
    hSigClone->SetMarkerSize(0.6);
    pad1->cd();
    hSigClone->Draw("E1 same");

    // 4) legend
    TLegend *leg = nullptr;
    {
        TIter next(pad1->GetListOfPrimitives());
        TObject *obj;
        while ((obj = next()))
            if (obj->InheritsFrom(TLegend::Class())) { leg = (TLegend*)obj; break; }
    }
    if (!leg) {
        leg = new TLegend(0.4, 0.2, 0.6, 0.4, NULL, "brNDC");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->Draw();
    }
    leg->AddEntry(hSigClone, signalLabel, "lep");

    pad1->Modified();
    cRatio->cd();
    cRatio->Modified();
    cRatio->Update();

    // 5) save
    gSystem->mkdir(outDir, kTRUE);
    TString stem = OutStem(ratioMacro);
    cRatio->SaveAs(Form("%s/%s__withSignal.pdf", outDir, stem.Data()));
    cRatio->SaveAs(Form("%s/%s__withSignal.C",   outDir, stem.Data()));

    std::cout << "Done: " << signalLabel << " overlaid on " << stem << std::endl;
}


void CombineHistos()
{
    //MET_trg_eff("../output/Gluino2000_Run2_METtrgEff_V11p15_Eta2p4.root", false);
    //MET_trg_eff("../output/Gluino2000_Run2_METtrgEff_AOD_V11p15_Eta2p4.root", true);
    //MET_trg_eff("METanalysis_Eta2p4_EffTrg", "../output/Gluino_V19/Gluino_Run3_MET_madgraph_2000_V19p1.root", false);
    //MET_trg_eff("METanalysis_Eta2p4_EffTrg", "../output/Gluino_V19/Gluino_Run3_MET_pythia_2000_V19p0.root", false);
    //PFMET_Cut(false);
    //TrigEff_AODvsMiniAOD();

    //MET_trg_eff("CalibPseudoMET", "MET_trg_eff_MC_vs_data", "../output/MuonEG_V17/MuonEG2024_V17p3.root", "../output/TTbar2024_V15/TTbar2024_V15p5.root");
    //MET_trg_eff("METanalysis_Eta2p4_EffTrg", "../output/Gluino_V13/Gluino2000_Run2_MET_V13p1.root", false);
    //Comp_muonEG("../output/MuonEG_V17/MuonEG2024_V17p3.root", "../output/TTbar2024_V15/TTbar2024_V15p5.root");
    // MET_trg_eff("CalibPseudoMET_MuWay", "CalibPseudoMET_MuWay", "WMuNu_PseudoMET", "../output/Mu2024_V18/Mu2024_V18.root", "../output/Wjets2024_V14/WjetMuNu2024_V14p8_weighted.root");
    // MET_trg_eff("CalibPseudoMET_MuWay", "CalibPseudoMET_MuWay_isRescaled", "WMuNu_PseudoMET_rescaled", "../output/Mu2024_V18/Mu2024_V18.root", "../output/Wjets2024_V14/WjetMuNu2024_V14p8_weighted.root");
    // MET_trg_eff("CalibPseudoMET", "CalibPseudoMET_isRescaled", "TTbar_PseudoMET_rescaled", "../output/MuonEG_V17/MuonEG2024_V17p4.root", "../output/TTbar2024_V15/TTbar2024_V15p8_weighted.root");
    // MET_trg_eff("CalibPseudoMET_MuWay", "CalibPseudoMET_MuWay_OTHERisRescaled", "TriggerEff_Mu2024_Wjets_PseudoMET_OTHERrescaled", "../output/Mu2024_V18/Mu2024_V18.root", "../output/Wjets2024_V14/WjetMuNu2024_V14p10_weighted.root");
    

    //Old_vs_New_fits("../output/JetMET2024_V12/JetMET2024_V12p24.root");
    // BKGdependency_preliminary("../output/TTbar2024_V15/TTbar2024_V15p5_weighted.root", "TTbar");
    // BKGdependency_preliminary("../output/Wjets2024_V14/Wjets2024_V14p6_weighted.root", "Wjets");
    // BKGdependency_preliminary("../output/QCD2024_V16/QCD2024_mu_V16p1_weighted.root", "QCD");
    //BKGdependency();

    //GluinoP_mass(true, true);
    //GluinoP_mass(false, true);
    CompareMaping();
    //Ihand1oP_fits_bkg();
    //etaSlices_bkg();

    //ShowSignalEfficiency();

    //FpixelPlot();
    //PseudoMET_vs_PFMET();
    //DrawPseudoMET(true); DrawPseudoMET(false);
    //nPV();
    //nHSCP();

    //ExtractSF ("TriggerEff_Mu2024_WMuNu_PseudoMETrescaled", "CalibPseudoMET_MuWay", "CalibPseudoMET_MuWay_isRescaled", "../output/Mu2024_V18/Mu2024_V18.root", "../output/Wjets2024_V14/WjetMuNu2024_V14p8_weighted.root");
    //ExtractSF ("TriggerEff_Mu2024_WMuNu_PseudoMET", "CalibPseudoMET_MuWay", "CalibPseudoMET_MuWay", "../output/Mu2024_V18/Mu2024_V18.root", "../output/Wjets2024_V14/WjetMuNu2024_V14p8_weighted.root");
    //ExtractSF ("TriggerEff_MuonEG2024_TTbar_PseudoMETrescaled", "CalibPseudoMET", "CalibPseudoMET_isRescaled", "../output/MuonEG_V17/MuonEG2024_V17p4.root", "../output/TTbar2024_V15/TTbar2024_V15p8.root");
    //PostTriggerPseudoMET();
    //PlotSF();

    //MyClusters();


    // Run2_vs_Run3_gluino__TriggerEff();
    // ComparePseudoMET(true); ComparePseudoMET(false);
    // YieldAfterSF(true, true, false, false); YieldAfterSF(false, true, false, false); // HLT Mu + HSCP selections; rescaled and not rescaled
    // YieldAfterSF(true, false, true, false); YieldAfterSF(false, false, true, false); // HLT MET + HSCP selections; rescaled and not rescaled
    // YieldAfterSF(true, false, false, true); YieldAfterSF(false, false, false, true); // HLT MuonEG + HSCP selections; rescaled and not rescaled
    // TableGluino(true); TableGluino(false);

    //LangausFitOnIh();

    //FpixSlices();
    //Ihand1oP_fits();
    //DefineIhCut();

    //SignalAcceptance_EtaSlice();
    //SignalAcceptance_EtaSlice(true);
    //Corr_Ih_1oP("Eta2p4"); Corr_Ih_1oP("Eta1"); Corr_Ih_1oP("Eta1_2p4");

    //Acceptance_EtaSlice_DataMC();


    //------------------------------------------------------------------
    // Trigger efficiencies
    //------------------------------------------------------------------

    // TriggerEffCalib__Signal("TriggerEffCalib__Signal",
    //                         "../output/Gluino_V19/Gluino_Run3_MET_madgraph_2000_V19p8.root", 
    //                         "Gluino2000");

    //TriggerEffCalib__Signal__2D("TriggerEffCalib__Signal", "../output/Gluino_V19/Gluino_Run3_MET_madgraph_2000_V19p4.root", "Gluino2000");

    // TriggerEffCalib("TriggerEffCalib", "TriggerEffCalib",
    //                 "../output/Mu2024_V18/Mu2024_V18p1.root", "../output/Wjets2024_V14/WjetMuNu2024_V14p12.root",
    //                 "TriggerEffCalib");

    // TriggerEffCalib__2D("TriggerEffCalib", "TriggerEffCalib",
    //                     "../output/Mu2024_V18/Mu2024_V18p1.root", "../output/Wjets2024_V14/WjetMuNu2024_V14p12.root",
    //                     "TriggerEffCalib");

    //SignalEffVsMass();
    // DisplayTriggerEff("../output/Gluino_V19/Gluino_Run3_MET_madgraph_2000_V19p8.root", true);
    // DisplayTriggerEff("../output/Mu2024_V18/Mu2024_V18p1.root", false);

    // ExtractSF("../outputDisplay/TriggEff/SF_PseudoMET.txt",
    //         "../outputDisplay/TriggEff/SF_PseudoMET_tex.txt",
    //         "TriggerEffCalib_PseudoCaloMET",
    //         "TriggerEffCalib_if___orMETtrg___PseudoCaloMET",
    //         "TriggerEffCalib_PseudoCaloMET",
    //         "TriggerEffCalib_if___orMETtrg___PseudoCaloMET",
    //         "../output/Mu2024_V18/Mu2024_V18p1.root",
    //         "../output/Wjets2024_V14/WjetMuNu2024_V14p12.root");

    // ExtractSF("../outputDisplay/TriggEff/SF_PUppiMET.txt",
    //         "../outputDisplay/TriggEff/SF_PUppiMET_tex.txt",
    //         "TriggerEffCalib_PUppiMET",
    //         "TriggerEffCalib_if___orMETtrg___PUppiMET",
    //         "TriggerEffCalib_PUppiMET",
    //         "TriggerEffCalib_if___orMETtrg___PUppiMET",
    //         "../output/Mu2024_V18/Mu2024_V18p1.root",
    //         "../output/Wjets2024_V14/WjetMuNu2024_V14p12.root");

    //FpixelInSignalAndData();

    //PairTypeStages_SingleMass(); PFType_ProportionAndTrigEff(); TriggerEfficiency_ByMass();

    // PairTypeStages_SingleMass__Notrigger();
    // PairTypeStages_SingleMass__Trigger();
    // TriggerEfficiency_ByMass();
    //TriggerEfficiency_VsBeta();

    //CompareKinematics();

    //ProfileVsRunNumber();

    //CompareIhWithCDF();

    // PUppiVSPF("../output/Gluino_V19/Gluino_Run3_MET_madgraph_2000_V19p9_weighted.root", false);
    //PUppiVSPF("../output/Mu2024_V18/Mu2024_V18p1.root", true);


    // Nm1Eff(false,
    //        "../output/TTbar2024_V15/TTbar2024_V15p9_weighted.root",
    //        "../output/TTbar2024_V15/TTbarSemiLep2024_V22p0_weighted.root",
    //        "../output/Wjets2024_V14/WjetMuNu2024_V14p13_weighted.root",
    //        "../output/QCD2024_V16/QCD2024_mu_V16p2_weighted.root",
    //        "../output/JetMET2024_V12/JetMET2024_V12p32.root",
    //        "../output/Gluino_V19/Gluino_Run3_MET_madgraph_2000_V19p8_weighted.root");

    // Nm1EffCutZeroed(false,
    //     "../output/TTbar2024_V15/TTbar2024_V15p9_weighted.root",
    //     "../output/Wjets2024_V14/WjetMuNu2024_V14p13_weighted.root",
    //     "../output/QCD2024_V16/QCD2024_mu_V16p2_weighted.root",
    //     "../output/JetMET2024_V12/JetMET2024_V12p32.root",
    //     "../output/Gluino_V19/Gluino_Run3_MET_madgraph_2000_V19p8_weighted.root"
    // );

    // NoselEff(false,
    //           "../output/TTbar2024_V15/TTbar2024_V15p9_weighted.root",
    //           "../output/TTbar2024_V15/TTbarSemiLep2024_V22p0_weighted.root",
    //           "../output/Wjets2024_V14/WjetMuNu2024_V14p13_weighted.root",
    //           "../output/QCD2024_V16/QCD2024_mu_V16p2_weighted.root",
    //           "../output/JetMET2024_V12/JetMET2024_V12p32.root",
    //           "../output/Gluino_V19/Gluino_Run3_MET_madgraph_2000_V19p8_weighted.root");

    // Cutflows("../output/QCD2024_V16/QCD2024_mu_V16p2_weighted.root",
    //          "../output/TTbar2024_V15/TTbar2024_V15p9_weighted.root",
    //          "../output/TTbar2024_V15/TTbarSemiLep2024_V22p0_weighted.root",
    //          "../output/Wjets2024_V14/WjetMuNu2024_V14p13_weighted.root",
    //          "../output/JetMET2024_V12/JetMET2024_V12p32.root", 
    //          "../output/Gluino_V19/Gluino_Run3_MET_madgraph_2000_V19p8_weighted.root");


    //PlotPthatQCD();
    //SignalMassResolution();

    //CompareSelEff("/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/Gluino_V19/", "SelEff");
    //PlotMassSpectraPerEta("/safe/ui3_1/cms/gcoulon/CMSSW_15_0_13_patch1/src/TupleAnalysis/output/Gluino_V19/", "SelEff", false, true, 1);

    FitVRRatios("../output/JetMET2024_V12/JetMET2024_V12p35_rebinEta4_rebinIh4_rebinP2_EtaReweighting_Eta1_OldFit_IhC.root", "eta1", true, true);
    FitVRRatios("../output/JetMET2024_V12/JetMET2024_V12p35_rebinEta4_rebinIh4_rebinP2_EtaReweighting_Eta1_2p4_OldFit_IhC.root", "eta1_2p4", true, true);
    FitVRRatios("../output/JetMET2024_V12/JetMET2024_V12p35_rebinEta4_rebinIh4_rebinP2_EtaReweighting_Eta2p4_OldFit_IhC.root", "eta2p4", true, true);

    return;
}