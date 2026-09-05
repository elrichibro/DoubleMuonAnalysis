#ifndef ANALYSISTOOLS_H
#define ANALYSISTOOLS_H

#include "TH1D.h"

struct  MC_Template_RooF{
    int pt_bin_idx;
    int eta_bin_idx;
    double eta_min, eta_max;
    double pt_min, pt_max;
    
    TH1D* h_pass{nullptr};
    TH1D* h_fail{nullptr};
};





#endif