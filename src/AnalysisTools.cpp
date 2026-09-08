#include "AnalysisTools.h"

#include <map>
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


#include <TCanvas.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <iostream>
#include <vector>

// ------------------------------------------------------------------------------------------------------------------------------------
// MC Template Maker
// ------------------------------------------------------------------------------------------------------------------------------------

std::vector<ROOT::RDF::RResultPtr<TH3D>> TemplateMaker(ROOT::RDF::RNode node, const config_struct& cfg, std::vector<ROOT::RDF::RResultPtr<TH2D>>& entry_map) {
    std::vector<ROOT::RDF::RResultPtr<TH3D>> container;
    
    ROOT::RDF::RNode node_hist = node;

    std::vector<float> pt_bins = cfg.analysis.pt_bins;
    std::vector<float> eta_bins = cfg.analysis.eta_bins;
    std::vector<float> mll_bins(cfg.analysis.mll_bins);

    float step = (120.0 - 60.0) / (cfg.analysis.mll_bins - 1.0);
    for (int i = 0; i < cfg.analysis.mll_bins; i++) {
        mll_bins[i] = 60.0 + (i * step);
    }

    std::string sample = cfg.general.dataset; 

    std::string name_pass = sample + "_h3_pass";
    std::string name_fail = sample + "_h3_fail";
    std::string name_entries_pass = sample + "_h2_entries_pass";
    std::string name_entries_fail = sample + "_h2_entries_fail";

    std::string title = "3D Histogram_" + sample;
    std::string title_entries = "2D Entries Histogram_" + sample;

    ROOT::RDF::TH3DModel model_pass(name_pass.c_str(), title.c_str(), eta_bins.size() - 1, eta_bins.data(), pt_bins.size() - 1, pt_bins.data(),
     mll_bins.size() - 1, mll_bins.data() );
    
    ROOT::RDF::TH3DModel model_fail(name_fail.c_str(), title.c_str(), eta_bins.size() - 1, eta_bins.data(), pt_bins.size() - 1, pt_bins.data(),
     mll_bins.size() - 1, mll_bins.data() );

    ROOT::RDF::TH2DModel model_entries_pass(name_entries_pass.c_str(), title_entries.c_str(), eta_bins.size() - 1, eta_bins.data(), pt_bins.size() - 1, pt_bins.data());
    ROOT::RDF::TH2DModel model_entries_fail(name_entries_fail.c_str(), title_entries.c_str(), eta_bins.size() - 1, eta_bins.data(), pt_bins.size() - 1, pt_bins.data());
    
    node_hist = node_hist
        .Define(sample + "_Probe_Pt_Pass", sample + "_Probe_Pt[" + sample + "_Mask_Pass]")
        .Define(sample + "_Probe_Eta_Pass", sample + "_Probe_Eta[" + sample + "_Mask_Pass]")
        .Define(sample + "_Mll_Pass", sample + "_Mll[" + sample + "_Mask_Pass]")
        .Define(sample + "_Probe_Pt_Fail", sample + "_Probe_Pt[!" + sample + "_Mask_Pass]")
        .Define(sample + "_Probe_Eta_Fail", sample + "_Probe_Eta[!" + sample + "_Mask_Pass]")
        .Define(sample + "_Mll_Fail", sample + "_Mll[!" + sample + "_Mask_Pass]");


    auto h3_pass = node_hist.Histo3D(model_pass, sample + "_Probe_Eta_Pass", sample + "_Probe_Pt_Pass", sample + "_Mll_Pass");
    auto h3_fail = node_hist.Histo3D(model_fail, sample + "_Probe_Eta_Fail", sample + "_Probe_Pt_Fail", sample + "_Mll_Fail");
    
    auto h2_pass = node_hist.Histo2D(model_entries_pass, sample + "_Probe_Eta_Pass", sample + "_Probe_Pt_Pass");
    auto h2_fail = node_hist.Histo2D(model_entries_fail, sample + "_Probe_Eta_Fail", sample + "_Probe_Pt_Fail");

    h2_pass->SetOption("COLZ TEXT");
    h2_fail->SetOption("COLZ TEXT");

    entry_map.push_back(h2_pass);
    entry_map.push_back(h2_fail);

    container.push_back(h3_pass);
    container.push_back(h3_fail);

    return container;
}

// ------------------------------------------------------------------------------------------------------------------------------------
// Template Loader
// ------------------------------------------------------------------------------------------------------------------------------------

int LoadTemplate(const config_struct& cfg, std::vector<Template_RooF>& container) {
    std::unique_ptr<TFile> file(TFile::Open(cfg.analysis.o_template_file_data.c_str(), "READ"));
    
    if (!file || file->IsZombie()) {
        std::cout << "ERROR: Cannot open the template output file: " << cfg.analysis.o_template_file_data << std::endl;
        return 1;
    }
    TDirectory* bin_dir = file->GetDirectory(cfg.analysis.bins_settup.c_str());
    
    if (!bin_dir) {
        std::cout << "ERROR: read operation fails, exiting." << std::endl;
        return 1;
    }

    auto h3_MC_pass = bin_dir->Get<TH3D>("MC_h3_pass");
    auto h3_MC_fail = bin_dir->Get<TH3D>("MC_h3_fail");
    auto h3_DATA_pass = bin_dir->Get<TH3D>("DATA_h3_pass");
    auto h3_DATA_fail = bin_dir->Get<TH3D>("DATA_h3_fail");

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
                bin_eta,
                bin_pt, 
                
                eta_bins[j], 
                eta_bins[j + 1],

                pt_bins[i], 
                pt_bins[i + 1], 
            
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
// Checker Pre-Fit Data
// ------------------------------------------------------------------------------------------------------------------------------------

void CheckPlotsTemplate(const std::vector<Template_RooF>& analysis_struct, const config_struct& cfg) {

    for (const auto& item : analysis_struct) {

        double eta_low  = cfg.analysis.eta_bins.at(item.eta_bin_idx - 1);
        double eta_high = cfg.analysis.eta_bins.at(item.eta_bin_idx);
        
        double pt_low   = cfg.analysis.pt_bins.at(item.pt_bin_idx - 1);
        double pt_high  = cfg.analysis.pt_bins.at(item.pt_bin_idx);

        std::string title = "( " + std::to_string(item.eta_bin_idx) + ", " + std::to_string(item.pt_bin_idx) +  ")";

        TCanvas* c_prefit = new TCanvas("c_prefit", title.c_str(), 1200, 800);
        c_prefit->Divide(2, 2);

        c_prefit->cd(1); 
        if (item.h_DATA_pass) {
            item.h_DATA_pass->Draw("HIST E");
        }
        c_prefit->cd(2); 
        if (item.h_DATA_fail) {
            item.h_DATA_fail->Draw("HIST E");
        }
        c_prefit->cd(3); 
        if (item.h_MC_pass) {
            item.h_MC_pass->Draw("HIST E");
        }
        c_prefit->cd(4); 
        if (item.h_MC_fail) {
            item.h_MC_fail->Draw("HIST E");
        }

        c_prefit->Update();

        while (gROOT->GetListOfCanvases()->FindObject("c_prefit")) {
            gSystem->ProcessEvents();
            gSystem->Sleep(50);
        }
    }
}

// ------------------------------------------------------------------------------------------------------------------------------------
// Fitter Binned data
// ------------------------------------------------------------------------------------------------------------------------------------

int Eff_BinnedFit(std::vector<Template_RooF>& analysis_struct, const analysis_params& params, std::vector<FitResult>& results, const int verb) {
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
        
        if (verb) {
            std::cout << "Events in Pass DATA: " << data_pass_entries <<std::endl;
            std::cout << "Events in Fail DATA: " << data_fail_entries <<std::endl;
            std::cout << "Events in Pass MC: " << mc_pass_entries <<std::endl;
            std::cout << "Events in Fail MC: " << mc_fail_entries <<std::endl;
        }  

        // SIGNAL
        RooRealVar efficiency("efficiency", "Efficiency", params.efficiency.at(0), params.efficiency.at(1), params.efficiency.at(2));
        RooRealVar n_sig_tot("n_sig_tot", "Number of signal events", total_entries * 0.9, 0.0, 1e8);
        
        // Gaussian convolution parameters
        RooRealVar mu("mu", "Mean gaussian", params.mu.at(0), params.mu.at(1), params.mu.at(2));
        RooRealVar sigma("sigma", "Sigma gaussina", params.sigma.at(0), params.sigma.at(1), params.sigma.at(2));
        RooGaussian gauss("gauss", "Smearing", mll, mu, sigma);

        RooFFTConvPdf conv_pass_sig("conv_pass_sig", "Conv model Pass signal + gauss", mll, hist_mc_pass_pdf, gauss);
        RooFFTConvPdf conv_fail_sig("conv_fail_sig", "Conv model Fail signal + gauss", mll, hist_mc_fail_pdf, gauss);

        // N events
        RooFormulaVar n_sig_pass("n_sig_pass", "efficiency * n_sig_tot", RooArgList(efficiency, n_sig_tot));
        RooFormulaVar n_sig_fail("n_sig_fail", "(1.0 - efficiency) * n_sig_tot", RooArgList(efficiency, n_sig_tot));

        // BACKGROUND
        RooRealVar lambda_pass("lambda_pass", "Decay freq Pass Bkg", params.lambda_pass.at(0), params.lambda_pass.at(1), params.lambda_pass.at(2));
        RooRealVar lambda_fail("lambda_fail", "Decay freq Fail Bkg", params.lambda_fail.at(0), params.lambda_fail.at(1), params.lambda_fail.at(2));
        
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

        /*
        Check reset
        mu.setVal(0.0);
        sigma.setVal(1.0);
        lambda_pass.setVal(-0.01);
        lambda_fail.setVal(-0.01);
        efficiency.setVal(0.90);
        */


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
            it.eta_bin_idx,
            it.pt_bin_idx,
            
            efficiency.getVal(),
            efficiency.getError(),
            
            n_sig_tot.getVal(),
            n_sig_tot.getError(),

            mu.getVal(),
            mu.getError(),
            
            sigma.getVal(),
            sigma.getError(),
                
            lambda_pass.getVal(),
            lambda_pass.getError(),

            lambda_fail.getVal(),
            lambda_fail.getError(),
            
            raw_fitRes ? raw_fitRes->status() : -1
        });
    }
    
    std::cout << "Succeded fits: " << successful_fits << std::endl;
    std::cout << "Failed fits: " << failed_fits << std::endl;

    std::cout << "Fit procedure finished." << std::endl;
    
    return 0;
}

// ------------------------------------------------------------------------------------------------------------------------------------
// Data Saver Post-Fit
// ------------------------------------------------------------------------------------------------------------------------------------

int SaveFitPlots(const std::vector<FitResult>& results, const config_struct& cfg, const std::vector<std::string>& vals) {
    const auto& eta_bins = cfg.analysis.eta_bins;
    const auto& pt_bins  = cfg.analysis.pt_bins;

    TFile o_template_file(cfg.analysis.o_template_file_data.c_str(), "UPDATE");
    if (o_template_file.IsZombie()) {
        std::cerr << "ERROR: Cannot open output file: " << cfg.analysis.o_template_file_data << std::endl;
        return 1;
    }

    std::string dir_name = "Plots_" + cfg.analysis.bins_settup;

    TDirectory* bin_dir = o_template_file.GetDirectory(dir_name.c_str());
    if (!bin_dir) {
        bin_dir = o_template_file.mkdir(dir_name.c_str());
    }

    bin_dir->cd();

    for (const auto& var_names : vals) {
        std::string histo_name  = "h2_" + var_names;
        std::string histo_title = var_names + " map;#eta;p_{T} [GeV]";

        TH2D h2_map(histo_name.c_str(), histo_title.c_str(), eta_bins.size() - 1, eta_bins.data(), pt_bins.size() - 1,  pt_bins.data());

        for (const auto& res : results) {
            double value = 0.0;
            double error = 0.0;
            
            if (var_names == "efficiency") {
                value = res.efficiency;
                error = res.efficiency_err;
            } else if (var_names == "n_tot") {
                value = res.n_tot;
                error = res.n_tot_err;
            
            } else if (var_names == "fit_status") {
                value = res.fit_status;

            } else if (var_names == "mu") {
                value = res.mu;
                error = res.mu_err;
            } else if (var_names == "sigma") {
                value = res.sigma;
                error = res.sigma_err;
            
            } else if (var_names == "lambda_pass") {
                value = res.lambda_pass;
                error = res.lambda_pass_err;
            } else if (var_names == "lambda_fail") {
                value = res.lambda_fail;
                error = res.lambda_fail_err;
            
            }  else {
                std::cout << "ERROR: Invalid fitted variable name: " << var_names << ". Exiting..." << std::endl;
                return 1;
            }

            if (var_names != "fit_status") {
                if (res.fit_status == 0) {
                    h2_map.SetBinContent(res.eta_bin_idx, res.pt_bin_idx, value);
                    h2_map.SetBinError(res.eta_bin_idx, res.pt_bin_idx, error);
                } else {
                    h2_map.SetBinContent(res.eta_bin_idx, res.pt_bin_idx, -999);
                    h2_map.SetBinError(res.eta_bin_idx, res.pt_bin_idx, -999);
                }
            } else {
                h2_map.SetBinContent(res.eta_bin_idx, res.pt_bin_idx, value);
                h2_map.SetBinError(res.eta_bin_idx, res.pt_bin_idx, error);
            }
        }

        h2_map.SetOption("TEXTE COLZ");
        h2_map.Write(histo_name.c_str(), TObject::kOverwrite);
    }
    
    o_template_file.Close();
    
    return 0;
}