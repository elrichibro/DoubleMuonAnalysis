#ifndef ANALYSISTOOLS_H
#define ANALYSISTOOLS_H

#include "TH1D.h"
#include "Config.h"
#include <ROOT/RDataFrame.hxx>

#include <string>


struct Template_RooF{
    int pt_bin_idx;
    int eta_bin_idx;
    
    float eta_min;
    float eta_max;
    float pt_min;
    float pt_max;
    
    TH1D* h_MC_pass{nullptr};
    TH1D* h_MC_fail{nullptr};
    TH1D* h_DATA_pass{nullptr};
    TH1D* h_DATA_fail{nullptr};
};

struct FitResult {
    int pt_bin_idx;
    int eta_bin_idx;
    
    double efficiency;
    double efficiency_err;
    double n_tot;
    double n_tot_err;
    double mu_0;
    double mu_0_err;
    double sigma;
    double sigma_err;
    
    int fit_status;
};

// ------------------------------------------------------------------------------------------------------------------------------------
// Template Maker
// ------------------------------------------------------------------------------------------------------------------------------------
std::vector<ROOT::RDF::RResultPtr<TH3D>> TemplateMaker(ROOT::RDF::RNode node, const config_struct& cfg, std::vector<ROOT::RDF::RResultPtr<TH2D>>& entry_map);

// ------------------------------------------------------------------------------------------------------------------------------------
// Load Template
// ------------------------------------------------------------------------------------------------------------------------------------
int LoadTemplate(const config_struct& cfg, std::vector<Template_RooF>& container);

// ------------------------------------------------------------------------------------------------------------------------------------
// EffBinnedFit
// ------------------------------------------------------------------------------------------------------------------------------------

int Eff_BinnedFit(std::vector<Template_RooF>& container, const config_struct& cfg, std::vector<FitResult>& results);

#endif