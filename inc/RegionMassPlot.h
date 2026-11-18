#ifndef RegionMassPlot_h
#define RegionMassPlot_h


#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <string>
#include <vector>

// Not needed by this class: HSCPSelector.h gets GetMass and loadSF1D/2D through this include.
#include "MassTools.h"

class RegionMassPlot{

    public:
        RegionMassPlot(std::string suffix,int etabins,int ihbins,int pbins,int massbins,int fpixbins, float C_parameter);
        ~RegionMassPlot();
       
        //Methods
        void initHisto(int& etabins,int& ihbins,int& pbins,int& massbins,int& fpixbins, float& C_parameter);
        void fill(float eta, float nhits, float p, float pt, float pterr, float ih, float ias, float m, float npv, float fpix, float w);
        void addToList(TList* list);	
        // The 26 histograms, in a fixed order (not the projections that addToList adds).
        std::vector<TH1*> histograms() const;
        // Adds the histograms of another plot with the same binning, histogram by histogram.
        void addFrom(const RegionMassPlot& other);

	//Data members
	std::string suffix_;
       
	//Plots binning
        int np;
        float plow;
        float pup;
        int npt;
        float ptlow;
        float ptup;
        int nih;
        float ihlow;
        float ihup;
        int nias;
        float iaslow;
        float iasup;
        int neta;
        float etalow;
        float etaup;
        int nmass;
        float masslow;
        float massup;
        int nfpix;
        float fpixlow;
        float fpixup;
        
        //List of all histos 
        TH2F* ih_pt;
        TH2F* ias_pt;
        TH2F* ih_ias;
        TH2F* ih_fpix;
        TH2F* eta_fpix;
        TH2F* oP_fpix;
        TH2F* ih_nhits;
        TH2F* ias_nhits;
        TH2F* eta_pt;
        TH2F* eta_1oP;
        TH2F* eta_p;
        TH2F* eta_pterrOpt;
        TH2F* nhits_pt;
        TH2F* eta_nhits;
        TH2F* ih_eta;
        TH2F* ih_p;
        TH2F* ias_p;
        TH2F* pt_pterroverpt;
        TH2F* ias_eta;
        TH2F* mass_eta;
        TH2F* eta_npv;
        TH2F* p_npv;
        TH2F* ih_npv;
        TH2F* mass_p;
        TH2F* mass_ih;
	TH1F* mass;
};

#endif
