#include "AnalysisTools.h"

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <memory>

#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooFormulaVar.h>
#include <RooGaussian.h>
#include <RooFFTConvPdf.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooSimultaneous.h>
#include <RooFitResult.h>

#include <TH2D.h>

// ------------------------------------------------------------------------------------------------------------------------------------
// MC Template Maker
// ------------------------------------------------------------------------------------------------------------------------------------

ROOT::RDF::RResultPtr<TH3D> TemplateMaker(ROOT::RDF::RNode node, const config_struct& cfg, const bool mask) {
    ROOT::RDF::RNode node_hist = node;

    std::vector<float> pt_bins = cfg.analysis.pt_bins;
    std::vector<float> eta_bins = cfg.analysis.eta_bins;
    std::vector<float> mll_bins(cfg.analysis.mll_bins);

    float step = (120.0 - 60.0) / (cfg.analysis.mll_bins - 1.0);
    for (int i = 0; i < cfg.analysis.mll_bins; i++) {
        mll_bins[i] = 60.0 + (i * step);
    }

    std::string sample = cfg.general.dataset; 

    std::string name = (mask) ? "h3_pass" : "h3_fail";
    name = sample + "_" + name;

    std::string title = "3D Histogram_" + sample;

    ROOT::RDF::TH3DModel model(name.c_str(), title.c_str(), eta_bins.size() - 1, eta_bins.data(), pt_bins.size() - 1, pt_bins.data(),
     mll_bins.size() - 1, mll_bins.data() );
    
    if (mask) {
        node_hist = node_hist
            .Define(sample + "_Probe_Pt_Pass", sample + "_Probe_Pt[" + sample + "_Mask_Pass]")
            .Define(sample + "_Probe_Eta_Pass", sample + "_Probe_Eta[" + sample + "_Mask_Pass]")
            .Define(sample + "_Mll_Pass", sample + "_Mll[" + sample + "_Mask_Pass]");

        auto h3 = node_hist.Histo3D(model, sample + "_Probe_Eta_Pass", sample + "_Probe_Pt_Pass", sample + "_Mll_Pass");
        return h3;
    } else {
        node_hist = node_hist
            .Define(sample + "_Probe_Pt_Fail", sample + "_Probe_Pt[!" + sample + "_Mask_Pass]")
            .Define(sample + "_Probe_Eta_Fail", sample + "_Probe_Eta[!" + sample + "_Mask_Pass]")
            .Define(sample + "_Mll_Fail", sample + "_Mll[!" + sample + "_Mask_Pass]");

        auto h3 = node_hist.Histo3D(model, sample + "_Probe_Eta_Fail", sample + "_Probe_Pt_Fail", sample + "_Mll_Fail");
        return h3;
    }
}


int LoadTemplate(const config_struct& cfg, std::vector<Template_RooF>& container) {
    std::unique_ptr<TFile> file(TFile::Open(cfg.io.o_file_template.c_str(), "READ"));
    
    if (!file || file->IsZombie()) {
        std::cout << "ERROR: Cannot open the template output file: " << cfg.io.o_file_template << std::endl;
        return 1;
    }

    auto h3_MC_pass = file->Get<TH3D>("MC_h3_pass");
    auto h3_MC_fail = file->Get<TH3D>("MC_h3_fail");
    auto h3_DATA_pass = file->Get<TH3D>("DATA_h3_pass");
    auto h3_DATA_fail = file->Get<TH3D>("DATA_h3_fail");

    if ((!h3_MC_pass) && (!h3_MC_fail) && (!h3_DATA_pass) && (!h3_DATA_fail)) {
        std::cout << "ERROR: null histogram pointer, exiting" << std::endl;
        return 1;
    }

    std::vector<float> pt_bins = cfg.analysis.pt_bins;
    std::vector<float> eta_bins = cfg.analysis.eta_bins;
    
    const int n_pt_bins = pt_bins.size() - 1;
    const int n_eta_bins = eta_bins.size() - 1;

    const int tot_bins = n_pt_bins * n_eta_bins;
    
    container.clear();
    container.reserve(tot_bins);

    std::cout << "Starting load..." << std::endl;

    for (int i = 0; i < n_pt_bins; i++) {
        int bin_pt = (i + 1);

        for (int j = 0; j < n_eta_bins; j++) {
            int bin_eta = (j + 1);
    
            std::string name_mc_pass = "h_mll_mc_pass_eta" + std::to_string(bin_eta) + "_pt" + std::to_string(bin_pt);
            std::string name_mc_fail = "h_mll_mc_fail_eta" + std::to_string(bin_eta) + "_pt" + std::to_string(bin_pt);
            std::string name_data_pass = "h_mll_data_pass_eta" + std::to_string(bin_eta) + "_pt" + std::to_string(bin_pt);
            std::string name_data_fail = "h_mll_data_fail_eta" + std::to_string(bin_eta) + "_pt" + std::to_string(bin_pt);
            
            TH1D* h1_mc_pass = h3_MC_pass->ProjectionZ(name_mc_pass.c_str(), bin_eta, bin_eta, bin_pt, bin_pt);
            TH1D* h1_mc_fail = h3_MC_fail->ProjectionZ(name_mc_fail.c_str(), bin_eta, bin_eta, bin_pt, bin_pt);
            TH1D* h1_data_pass = h3_DATA_pass->ProjectionZ(name_data_pass.c_str(), bin_eta, bin_eta, bin_pt, bin_pt);
            TH1D* h1_data_fail = h3_DATA_fail->ProjectionZ(name_data_fail.c_str(), bin_eta, bin_eta, bin_pt, bin_pt);
        
            if (h1_mc_pass) {
                h1_mc_pass->SetDirectory(nullptr);
            }
            if (h1_mc_fail) {
                h1_mc_fail->SetDirectory(nullptr);
            }
            if (h1_data_pass) {
                h1_data_pass->SetDirectory(nullptr);
            }
            if (h1_data_fail) {
                h1_data_fail->SetDirectory(nullptr);
            }
            
            container.emplace_back(Template_RooF{
                bin_pt, 
                bin_eta, 
                pt_bins[i], 
                pt_bins[i + 1], 
                eta_bins[j], 
                eta_bins[j + 1],
        
                h1_mc_pass,
                h1_mc_fail,
                h1_data_pass,
                h1_data_fail
            });
        }
    }
    std::cout << "Load finished, exiting..." << std::endl;

    file->Close();
    
    return 0;
}

// ------------------------------------------------------------------------------------------------------------------------------------
// ROOFIT
// ------------------------------------------------------------------------------------------------------------------------------------

int Eff_BinnedFit(std::vector<Template_RooF>& analysis_struct, const config_struct& cfg, std::vector<FitResult>& results) {
    results.reserve(analysis_struct.size());
    
    RooRealVar mll("mll", "m_{#mu+#mu-}", 60.0, 120.0, "GeV");
    mll.setBins(10000, "cache");

    RooCategory sample("sample", "TagAndProbe categories");
    sample.defineType("Pass");
    sample.defineType("Fail");

    int successful_fits = 0;
    int failed_fits = 0;

    for (const auto& it : analysis_struct) {
        std::cout << "Fitting -> Eta bin: " << it.eta_bin_idx << " , Pt bin: " << it.pt_bin_idx << std::endl;

        RooDataHist hist_mc_pass("hist_mc_pass", "Pass data histogram", mll, it.h_MC_pass);
        RooDataHist hist_mc_fail("hist_mc_fail", "Fail data histogram", mll, it.h_MC_fail);

        RooHistPdf hist_mc_pass_pdf("hist_mc_pass_pdf", "Pass data pdf", mll, hist_mc_pass);
        RooHistPdf hist_mc_fail_pdf("hist_mc_fail_pdf", "Fail data pdf", mll, hist_mc_fail);

        std::map<std::string, TH1*> map_hist_data;
        map_hist_data["Pass"] = it.h_DATA_pass;
        map_hist_data["Fail"] = it.h_DATA_fail;

        RooDataHist hist_sig_data("hist_sig_data", "Signal data", mll, sample, map_hist_data);

        double data_pass_entries = it.h_DATA_pass->Integral();
        double data_fail_entries = it.h_DATA_fail->Integral();
        double total_entries = data_pass_entries + data_fail_entries;

        double mc_pass_entries = it.h_MC_pass->Integral();
        double mc_fail_entries = it.h_MC_fail->Integral();
        
        if (cfg.general.verbose) {
            std::cout << "Events in Pass DATA: " << data_pass_entries <<std::endl;
            std::cout << "Events in Fail DATA: " << data_fail_entries <<std::endl;
            std::cout << "Events in Pass MC: " << mc_pass_entries <<std::endl;
            std::cout << "Events in Fail MC: " << mc_fail_entries <<std::endl;
        }  

        // SIGNAL
        RooRealVar efficiency("efficiency", "Efficiency", 0.9, 0.5, 1.0);
        RooRealVar n_sig_tot("n_sig_tot", "Number of signal events", total_entries * 0.9, 0.0, 1e8);
        
        // Gaussian convolution parameters
        RooRealVar mu("mu", "Mean gaussian", 0.0, -5.0, 5.0);
        RooRealVar sigma("sigma", "Sigma gaussina", 1.0, 0.0001, 10.0);
        RooGaussian gauss("gauss", "Smearing", mll, mu, sigma);

        RooFFTConvPdf conv_pass_sig("conv_pass_sig", "Conv model Pass signal + gauss", mll, hist_mc_pass_pdf, gauss);
        RooFFTConvPdf conv_fail_sig("conv_fail_sig", "Conv model Fail signal + gauss", mll, hist_mc_fail_pdf, gauss);

        // N events
        RooFormulaVar n_sig_pass("n_sig_pass", "efficiency * n_sig_tot", RooArgList(efficiency, n_sig_tot));
        RooFormulaVar n_sig_fail("n_sig_fail", "(1.0 - efficiency) * n_sig_tot", RooArgList(efficiency, n_sig_tot));

        // BACKGROUND
        RooRealVar lambda_pass("lambda_pass", "Decay freq Pass Bkg", -0.01, -10.0, 0.0);
        RooRealVar lambda_fail("lambda_fail", "Decay freq Fail Bkg", -0.01, -10.0, 0.0);
        
        RooExponential bkg_pass("bkg_pass", "Bkg Pass", mll, lambda_pass);
        RooExponential bkg_fail("bkg_fail", "Bkg Fail", mll, lambda_fail);

        // N events
        RooRealVar n_bkg_pass("n_bkg_pass", "Number of Bkg Pass events", data_pass_entries * 0.1, 0.0, 1e7);
        RooRealVar n_bkg_fail("n_bkg_fail", "Number of Bkg Fail events", data_fail_entries * 0.5, 0.0, 1e7);


        // Global Model
        RooAddPdf model_pass("model_pass", "Global Pass model", RooArgList(conv_pass_sig, bkg_pass), RooArgList(n_sig_pass, n_bkg_pass));
        RooAddPdf model_fail("model_fail", "Global Fail model", RooArgList(conv_fail_sig, bkg_fail), RooArgList(n_sig_fail, n_bkg_fail));

        RooSimultaneous simPdf("simPdf", "Simultaneous Fit Pass/Fail stats", sample);
        simPdf.addPdf(model_pass, "Pass");
        simPdf.addPdf(model_fail, "Fail");


        std::unique_ptr<RooFitResult> raw_fitRes(simPdf.fitTo(
        
            hist_sig_data, 
            RooFit::Extended(true),// Extended LikelyHood -> fundamental
            RooFit::Save(true), 
            RooFit::PrintLevel(-1), 
            RooFit::PrintEvalErrors(-1),
            RooFit::Verbose(false)
        
        ));

        int fit_status = raw_fitRes ? raw_fitRes->status() : -1;
        if (fit_status == 0) {
            successful_fits++;
        } else {
            failed_fits++;
        }

        raw_fitRes->Print("v");

        results.emplace_back( FitResult {
            it.pt_bin_idx,
            it.eta_bin_idx,
            efficiency.getVal(),
            efficiency.getError(),
            n_sig_tot.getVal(),
            n_sig_tot.getError(),
            mu.getVal(),
            mu.getError(),
            sigma.getVal(),
            sigma.getError(),
            raw_fitRes ? raw_fitRes->status() : -1
        });
    }

    std::cout << "Fit procedure finished." << std::endl;
    return 0;
}