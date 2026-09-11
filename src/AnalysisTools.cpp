#include "AnalysisTools.h"

#include "TH1D.h"

#include <map>
#include <memory>

#include <RooHistPdf.h>
#include <RooFormulaVar.h>

#include <RooGaussian.h>
#include <RooExponential.h>
#include <RooAddPdf.h>

#include <RooFFTConvPdf.h>
#include <RooKeysPdf.h>

#include <RooFitResult.h>
#include <RooMinimizer.h>

#include <RooPlot.h>
#include <TPaveText.h>

#include <TCanvas.h>
#include <TSystem.h>

// ------------------------------------------------------------------------------------------------------------------------------------

int TemplateMaker(ROOT::RDF::RNode node, const config_struct& cfg, const int dataset) {   
    ROOT::RDF::RNode node_hist = node;

    std::vector<float> pt_bins = cfg.templ.pt_bins;
    std::vector<float> eta_bins = cfg.templ.eta_bins;
    std::vector<float> mll_bins(cfg.templ.mll_bins + 1);

    float step = (120.0 - 60.0) / cfg.templ.mll_bins;
    
    for (int i = 0; i <= cfg.templ.mll_bins; i++) {
        mll_bins[i] = 60.0 + (i * step);
    }

    std::string sample = "";
    
    if (dataset == 1) {
        sample = "DATA";
    } else if (dataset == 2) {
        sample = "MC";
    } else {
        std::cout << "ERROR: invalid dataset identifier, exiting..." << std::endl;
        return 1;
    }

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

    bool histo = false;
    bool snapshot = false;

    ROOT::RDF::RResultPtr<TH3D> h3_pass, h3_fail;
    ROOT::RDF::RResultPtr<TH2D> h2_pass, h2_fail;

    snapshot_type snap;

    if (cfg.templ.template_type.find("HISTO") != std::string::npos) {
        
        h3_pass = node_hist.Histo3D(model_pass, sample + "_Probe_Eta_Pass", sample + "_Probe_Pt_Pass", sample + "_Mll_Pass");
        h3_fail = node_hist.Histo3D(model_fail, sample + "_Probe_Eta_Fail", sample + "_Probe_Pt_Fail", sample + "_Mll_Fail");
        
        h2_pass = node_hist.Histo2D(model_entries_pass, sample + "_Probe_Eta_Pass", sample + "_Probe_Pt_Pass");
        h2_fail = node_hist.Histo2D(model_entries_fail, sample + "_Probe_Eta_Fail", sample + "_Probe_Pt_Fail");

        histo = true;
    }
    
    if (cfg.templ.template_type.find("DATA") != std::string::npos) {
        
        std::vector<std::string> names;

        node_hist = ApplyKinematicalBinDivision(node_hist, cfg, sample + "_Probe_Pt_Pass", sample + "_Probe_Eta_Pass", sample + "_Mll_Pass", 
            names, dataset, 1);
        
        std::cout << "First filter applied" << std::endl;

        node_hist = ApplyKinematicalBinDivision(node_hist, cfg, sample + "_Probe_Pt_Fail", sample + "_Probe_Eta_Fail", sample + "_Mll_Fail", 
            names, dataset, 2);

        std::cout << "Second filter applied" << std::endl;
        std::cout << "Initializing SnapShot" << std::endl;
        
        ROOT::RDF::RSnapshotOptions snapshot_opts;
        snapshot_opts.fMode = "UPDATE";
        snapshot_opts.fLazy = true;
        snapshot_opts.fOverwriteIfExists = true;


        std::string snaphot_name = sample + "_" + cfg.templ.bins_settup + "_Tree";
        snap = node_hist.Snapshot(snaphot_name, cfg.templ.o_template_file_data.c_str(), names, snapshot_opts);

        snapshot = true;
    }

    if (snapshot) {
        snap.GetValue(); 
    }
    
    if (histo) {
        std::cout << "Writing Histograms into " << cfg.templ.o_template_file_data << " file." << std::endl;

        TFile o_template_file(cfg.templ.o_template_file_data.c_str(), "UPDATE");

        if (o_template_file.IsZombie()) {
            std::cout << "ERROR: Cannot find output file: " << cfg.templ.o_template_file_data << std::endl;
            return 1;
        }

        TDirectory* bin_dir = o_template_file.GetDirectory(cfg.templ.bins_settup.c_str());
        if (!bin_dir) {
            bin_dir = o_template_file.mkdir(cfg.templ.bins_settup.c_str());
        }

        bin_dir->cd();    

        h2_pass->SetOption("COLZ TEXT");
        h2_fail->SetOption("COLZ TEXT");

        h3_pass->Write(nullptr, TObject::kOverwrite);
        h3_fail->Write(nullptr, TObject::kOverwrite);

        h2_pass->Write(nullptr, TObject::kOverwrite);
        h2_fail->Write(nullptr, TObject::kOverwrite);

        o_template_file.Write();
        o_template_file.Close();
    }

    std::cout << "Saving data in intermediate step as: " << cfg.templ.template_type << " in the file " << cfg.templ.o_template_file_data << std::endl;

    return 0;
}

// ------------------------------------------------------------------------------------------------------------------------------------
// Reading RVec into RooDataSet var
// ------------------------------------------------------------------------------------------------------------------------------------

RooDataSet* LoadRVecIntoDataset(TTree* tree, const std::string& branch_name, RooRealVar& mll, const std::string& data_name) {
    if (!tree) {
        std::cout << "ERROR: invalid tre number, exiting..." << std::endl;
        return nullptr;
    }

    RooDataSet* data = new RooDataSet(data_name.c_str(), "Unbinned mll", RooArgSet(mll));
 
    ROOT::VecOps::RVec<float>* column_vec = nullptr;
    tree->SetBranchAddress(branch_name.c_str(), &column_vec);
 
    Long64_t n_entries = tree->GetEntries();

    for (Long64_t i = 0; i < n_entries; i++) {
        tree->GetEntry(i);
 
        if (column_vec != nullptr) {
            for (float mass : *column_vec) {
                mll.setVal(mass);
                data->add(RooArgSet(mll));
            }
        }
    }
    tree->ResetBranchAddresses();

    return data;
}

/*
RooDataSet* LoadRVecIntoDataset(TTree* tree, const std::string& branch_name, RooRealVar& mll, const std::string& data_name) {
    if (!tree) return nullptr;

    TTree flat_tree("flat_tree", "Tempo tree");

    float mass_val = 0.0f;
    
    flat_tree.Branch(mll.GetName(), &mass_val, Form("%s/F", mll.GetName()));

    ROOT::VecOps::RVec<float>* column_vec = nullptr;
    tree->SetBranchAddress(branch_name.c_str(), &column_vec);

    Long64_t n_entries = tree->GetEntries();

    for (Long64_t i = 0; i < n_entries; i++) {
        tree->GetEntry(i);
        if (column_vec) {
            for (float mass : *column_vec) {
                mass_val = mass;
                flat_tree.Fill(); 
            }
        }
    }
    tree->ResetBranchAddresses();

    RooDataSet* data = new RooDataSet(
        data_name.c_str(), 
        "Unbinned mll", 
        RooArgSet(mll),
        RooFit::Import(flat_tree)
    );

    return data;
}
*/

// ------------------------------------------------------------------------------------------------------------------------------------
// Template Loader
// ------------------------------------------------------------------------------------------------------------------------------------

int LoadTemplate(const config_struct& cfg, std::vector<Template_RooF>& container) {
    RooRealVar mll("mll", "m_{#mu+#mu-}", 60.0, 120.0, "GeV");

    TH3D* h3_MC_pass = nullptr;
    TH3D* h3_MC_fail = nullptr;
    TH3D* h3_DATA_pass = nullptr;
    TH3D* h3_DATA_fail = nullptr;

    std::string name_data_fail = "";
    std::string name_data_pass = "";
    std::string name_mc_fail = "";
    std::string name_mc_pass = "";
  
    bool h_mc_pass = false;
    bool h_mc_fail = false;
    bool h_data_pass = false;
    bool h_data_fail = false;

    bool d_mc_pass = false;
    bool d_mc_fail = false;
    bool d_data_pass = false;
    bool d_data_fail = false;

    if ((cfg.analysis.sample_pass_mc == "histo")) {
        h_mc_pass = true;
    } else if((cfg.analysis.sample_pass_mc == "data")) {
        d_mc_pass = true;
    }

    if ((cfg.analysis.sample_pass_data == "histo")) {
        h_data_pass = true;
    } else if((cfg.analysis.sample_pass_data == "data")) {
        d_data_pass = true;
    }

    if ((cfg.analysis.sample_fail_mc == "histo")) {
        h_mc_fail = true;
    } else if((cfg.analysis.sample_fail_mc == "data")) {
        d_mc_fail = true;
    }

    if ((cfg.analysis.sample_fail_data == "histo")) {
        h_data_fail = true;
    } else if((cfg.analysis.sample_fail_data == "data")) {
        d_data_fail = true;
    }

    std::unique_ptr<TFile> file(TFile::Open(cfg.templ.o_template_file_data.c_str(), "READ"));

    if ((!file) || (file->IsZombie())) {
        std::cout << "ERROR: Cannot open the template output file: " << cfg.templ.o_template_file_data << std::endl;
        return 1;
    }

    TDirectory* bin_dir = nullptr;
    TTree* tree_data = nullptr;
    TTree* tree_mc = nullptr;
    
    if ((h_data_pass) || (h_data_fail) || (h_mc_fail) || (h_mc_pass)) {
        bin_dir = file->GetDirectory(cfg.templ.bins_settup.c_str());
        
        if (!bin_dir) {
            std::cout << "ERROR: read operation fails, exiting." << std::endl;
            return 1;
        }
    }

    if (d_data_pass || d_data_fail) {
        std::string tree_name = "DATA_" + cfg.templ.bins_settup + "_Tree";
        tree_data = file->Get<TTree>(tree_name.c_str());
    }

    if (d_mc_pass || d_mc_fail) {
        std::string tree_name = "MC_" + cfg.templ.bins_settup + "_Tree";
        tree_mc = file->Get<TTree>(tree_name.c_str());
    }

    // ----------------
    // Histogram option
    // ----------------

    if (h_mc_pass) {
        h3_MC_pass = bin_dir->Get<TH3D>("MC_h3_pass");
    }

    if (h_data_pass) {
        h3_DATA_pass = bin_dir->Get<TH3D>("DATA_h3_pass");
    }

    if (h_mc_fail) {
        h3_MC_fail = bin_dir->Get<TH3D>("MC_h3_fail");
    }

    if (h_data_fail) {
        h3_DATA_fail = bin_dir->Get<TH3D>("DATA_h3_fail");
    }

    // Getting bins
    std::vector<float> pt_bins = cfg.templ.pt_bins;
    std::vector<float> eta_bins = cfg.templ.eta_bins;
    
    const int n_pt_bins = pt_bins.size() - 1;
    const int n_eta_bins = eta_bins.size() - 1;

    const int tot_bins = n_pt_bins * n_eta_bins;
    
    // Reserve capacity
    container.clear();
    container.reserve(tot_bins);

    std::cout << "Starting load..." << std::endl;

        for (int i = 0; i < n_pt_bins; i++) {
            int bin_pt = (i + 1);

            for (int j = 0; j < n_eta_bins; j++) {
                int bin_eta = (j + 1);
                
                // New pointers
                TH1D* h1_mc_pass = nullptr;
                TH1D* h1_mc_fail = nullptr;
                TH1D* h1_data_pass = nullptr;
                TH1D* h1_data_fail = nullptr;

                RooDataSet* d1_mc_pass = nullptr;
                RooDataSet* d1_mc_fail = nullptr;
                RooDataSet* d1_data_pass = nullptr;
                RooDataSet* d1_data_fail = nullptr;

                // Histograms

                if (h_mc_pass) {
                    name_mc_pass = "h_mll_mc_pass_eta" + std::to_string(bin_eta) + "_pt" + std::to_string(bin_pt);
                    h1_mc_pass = h3_MC_pass->ProjectionZ(name_mc_pass.c_str(), bin_eta, bin_eta, bin_pt, bin_pt);
                    h1_mc_pass->SetDirectory(nullptr);
                }
                
                if (h_mc_fail) {
                    name_mc_fail = "h_mll_mc_fail_eta" + std::to_string(bin_eta) + "_pt" + std::to_string(bin_pt);
                    h1_mc_fail = h3_MC_fail->ProjectionZ(name_mc_fail.c_str(), bin_eta, bin_eta, bin_pt, bin_pt);
                    h1_mc_fail->SetDirectory(nullptr);
                }
                
                if (h_data_pass) {
                    name_data_pass = "h_mll_data_pass_eta" + std::to_string(bin_eta) + "_pt" + std::to_string(bin_pt);
                    h1_data_pass = h3_DATA_pass->ProjectionZ(name_data_pass.c_str(), bin_eta, bin_eta, bin_pt, bin_pt);
                    h1_data_pass->SetDirectory(nullptr);
                }

                if (h_data_fail) {
                    name_data_fail = "h_mll_data_fail_eta" + std::to_string(bin_eta) + "_pt" + std::to_string(bin_pt);
                    h1_data_fail = h3_DATA_fail->ProjectionZ(name_data_fail.c_str(), bin_eta, bin_eta, bin_pt, bin_pt);
                    h1_data_fail->SetDirectory(nullptr);
                }

                // Data

                if (d_mc_pass) {
                    std::string leaf_name = "MC_Mll_Eta" + std::to_string(bin_eta) + "_Pt" + std::to_string(bin_pt) + "_Pass";
                    d1_mc_pass = LoadRVecIntoDataset(tree_mc, leaf_name, mll, "d_" + leaf_name);
                }
                
                if (d_mc_fail) {
                    std::string leaf_name = "MC_Mll_Eta" + std::to_string(bin_eta) + "_Pt" + std::to_string(bin_pt) + "_Fail";
                    d1_mc_fail = LoadRVecIntoDataset(tree_mc, leaf_name, mll, "d_" + leaf_name);
                }
                
                if (d_data_pass) {
                    std::string leaf_name = "DATA_Mll_Eta" + std::to_string(bin_eta) + "_Pt" + std::to_string(bin_pt) + "_Pass";
                    d1_data_pass = LoadRVecIntoDataset(tree_data, leaf_name, mll, "d_" + leaf_name);
                }
                
                if (d_data_fail) {
                    std::string leaf_name = "DATA_Mll_Eta" + std::to_string(bin_eta) + "_Pt" + std::to_string(bin_pt) + "_Fail";
                    d1_data_fail = LoadRVecIntoDataset(tree_data, leaf_name, mll, "d_" + leaf_name);
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
                    h1_data_fail,

                    d1_mc_pass,
                    d1_mc_fail,
                    
                    d1_data_pass,
                    d1_data_fail
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
        
        double eta_low  = cfg.templ.eta_bins.at(item.eta_bin_idx - 1);
        double eta_high = cfg.templ.eta_bins.at(item.eta_bin_idx);
        
        double pt_low   = cfg.templ.pt_bins.at(item.pt_bin_idx - 1);
        double pt_high  = cfg.templ.pt_bins.at(item.pt_bin_idx);

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

int Eff_BinnedFit(std::vector<Template_RooF>& analysis_struct, const config_struct& cfg, std::vector<FitResult>& results, TFile* o_file) {
    
    results.reserve(analysis_struct.size());
    
    RooRealVar mll("mll", "m_{#mu+#mu-}", 60.0, 120.0, "GeV");
    mll.setBins(8192, "cache");// For FFT calculus

    RooCategory sample("sample", "TagAndProbe categories");
    sample.defineType("Pass");
    sample.defineType("Fail");

    int successful_fits = 0;
    int failed_fits = 0;

    for (const auto& it : analysis_struct) {
        std::cout << "Fitting -> Eta bin: " << it.eta_bin_idx << " , Pt bin: " << it.pt_bin_idx << std::endl;

        // --------------
        // MC PASS option
        // --------------

        std::unique_ptr<RooAbsPdf> mc_pass_pdf = nullptr;
        std::unique_ptr<RooDataHist> hist_mc_pass = nullptr;

        if (cfg.analysis.sample_pass_mc == "data") {
            mc_pass_pdf = std::make_unique<RooKeysPdf>("mc_pass_pdf", "MC Pass KeysPdf", mll, *it.d_MC_pass, RooKeysPdf::NoMirror, 1.5);

        } else if (cfg.analysis.sample_pass_mc == "histo") {

            hist_mc_pass = std::make_unique<RooDataHist>("hist_mc_pass", "Pass data histogram", mll, it.h_MC_pass);
            auto hist_pdf = std::make_unique<RooHistPdf>("hist_mc_pass_pdf", "Pass data pdf", mll, *hist_mc_pass);
            hist_pdf->setInterpolationOrder(2);
            mc_pass_pdf = std::move(hist_pdf);
        }

        // --------------
        // MC FAIL option
        // --------------

        std::unique_ptr<RooAbsPdf> mc_fail_pdf = nullptr;
        std::unique_ptr<RooDataHist> hist_mc_fail = nullptr;

        if (cfg.analysis.sample_fail_mc == "data") {
            mc_fail_pdf = std::make_unique<RooKeysPdf>("mc_fail_pdf", "MC Fail KeysPdf", mll, *it.d_MC_fail, RooKeysPdf::NoMirror, 1.5);

        } else if (cfg.analysis.sample_fail_mc == "histo") {

            hist_mc_fail = std::make_unique<RooDataHist>("hist_mc_fail", "Pass fail histogram", mll, it.h_MC_fail);
            auto hist_pdf = std::make_unique<RooHistPdf>("hist_mc_fail_pdf", "Pass fail pdf", mll, *hist_mc_fail);
            hist_pdf->setInterpolationOrder(2);
            mc_fail_pdf = std::move(hist_pdf);
        }

        // DATA container
        std::unique_ptr<RooAbsData> sig_data;

        if ((cfg.analysis.sample_pass_data == "data") && (cfg.analysis.sample_fail_data == "data")) {
    
            sig_data = std::make_unique<RooDataSet>("sig_data", "Signal data unbinned", RooArgSet(mll), 
            RooFit::Index(sample),
            RooFit::Import("Pass", *it.d_DATA_pass), 
            RooFit::Import("Fail", *it.d_DATA_fail));

        } else if ((cfg.analysis.sample_fail_data == "histo") && (cfg.analysis.sample_fail_data == "histo")){
        
            std::map<std::string, TH1*> map_hist_data;
            map_hist_data["Pass"] = it.h_DATA_pass;
            map_hist_data["Fail"] = it.h_DATA_fail;

            sig_data = std::make_unique<RooDataHist>("sig_data", "Signal data binned", mll, sample, map_hist_data);
        
        } else {
            std::cout << "ERROR: invalid combination in JSON file, exiting..." << std::endl;
            return 1;
        }


        bool is_data_unbinned = (cfg.analysis.sample_pass_data == "data");

        double data_pass_entries = is_data_unbinned ? it.d_DATA_pass->sumEntries() : it.h_DATA_pass->Integral();
        double data_fail_entries = is_data_unbinned ? it.d_DATA_fail->sumEntries() : it.h_DATA_fail->Integral();
        double total_entries = data_pass_entries + data_fail_entries;
        
        if (cfg.general.verbose) {
            std::cout << "Events in Pass DATA: " << data_pass_entries <<std::endl;
            std::cout << "Events in Fail DATA: " << data_fail_entries <<std::endl;
        }  


        // Efficiency
        RooRealVar efficiency("efficiency", "Efficiency", cfg.analysis.params.efficiency.at(0), cfg.analysis.params.efficiency.at(1),
         cfg.analysis.params.efficiency.at(2));
        // Number of signal events
        RooRealVar n_sig_tot("n_sig_tot", "Number of signal events", total_entries * 0.9, 0.0, 1e8);
        
        // Gaussian convolution parameters
        RooRealVar mu("mu", "Mean gaussian", cfg.analysis.params.mu.at(0), cfg.analysis.params.mu.at(1), cfg.analysis.params.mu.at(2));
        RooRealVar sigma("sigma", "Sigma gaussina", cfg.analysis.params.sigma.at(0), cfg.analysis.params.sigma.at(1), cfg.analysis.params.sigma.at(2));
        RooGaussian gauss("gauss", "Smearing", mll, mu, sigma);

        RooFFTConvPdf conv_pass_sig("conv_pass_sig", "Conv model Pass signal + gauss", mll, *mc_pass_pdf, gauss);
        RooFFTConvPdf conv_fail_sig("conv_fail_sig", "Conv model Fail signal + gauss", mll, *mc_fail_pdf, gauss);

        conv_pass_sig.setBufferFraction(0.2);
        conv_fail_sig.setBufferFraction(0.2);

        // N events
        RooFormulaVar n_sig_pass("n_sig_pass", "efficiency * n_sig_tot", RooArgList(efficiency, n_sig_tot));
        RooFormulaVar n_sig_fail("n_sig_fail", "(1.0 - efficiency) * n_sig_tot", RooArgList(efficiency, n_sig_tot));

        // BACKGROUND
        RooRealVar lambda_pass("lambda_pass", "Decay freq Pass Bkg", cfg.analysis.params.lambda_pass.at(0), cfg.analysis.params.lambda_pass.at(1), 
        cfg.analysis.params.lambda_pass.at(2));
        
        RooRealVar lambda_fail("lambda_fail", "Decay freq Fail Bkg", cfg.analysis.params.lambda_fail.at(0), cfg.analysis.params.lambda_fail.at(1), 
        cfg.analysis.params.lambda_fail.at(2));
        
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

        if (cfg.analysis.pre_fit) {
            efficiency.setVal(1.0);
            efficiency.setConstant(kTRUE);
            RooDataHist hist_data_pass("hist_data_pass", "Pass data histogram", mll, it.h_DATA_pass);

            std::unique_ptr<RooAbsReal> nll_pass(model_pass.createNLL(
                hist_data_pass,
                RooFit::Extended(kTRUE), 
                RooFit::NumCPU(1)
            ));

            RooMinimizer m_pass(*nll_pass);

            m_pass.setStrategy(2);
            m_pass.setEps(1e-2);
            int status_migrad = m_pass.migrad();
        
            if (status_migrad == 0 || status_migrad == 1) {
                m_pass.hesse();
            }
            //m_pass.minos();

            std::unique_ptr<RooFitResult> prefitRes(m_pass.save());

            double mu_fit = mu.getVal();
            double sigma_fit = sigma.getVal();
            double mu_err = mu.getError();
            double sigma_err = sigma.getError();

            mu.setConstant(kFALSE);
            sigma.setConstant(kFALSE);

            mu.setVal(mu_fit);
            sigma.setVal(sigma_fit);
        
            efficiency.setConstant(kFALSE);
            efficiency.setVal(cfg.analysis.params.efficiency.at(0));
        }

        std::unique_ptr<RooAbsReal> nll(simPdf.createNLL(
            *sig_data, 
            RooFit::Extended(kTRUE), 
            RooFit::NumCPU(1)
        )); 

        RooMinimizer m(*nll);

        m.setStrategy(2);
        m.setEps(1e-2);
        int status_migrad = m.migrad();
        
        if (status_migrad == 0 || status_migrad == 1) {
            m.hesse();
        }
        //m.minos();

        std::unique_ptr<RooFitResult> raw_fitRes(m.save());

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

        SaveBinFitCanvas(mll, sample, *sig_data, simPdf, results.back(), o_file, cfg);
    }
    
    std::cout << "Succeded fits: " << successful_fits << std::endl;
    std::cout << "Failed fits: " << failed_fits << std::endl;

    std::cout << "Fit procedure finished." << std::endl;
    
    return 0;
}

// ------------------------------------------------------------------------------------------------------------------------------------
// Data Saver Post-Fit
// ------------------------------------------------------------------------------------------------------------------------------------

int SaveMapFittedValues(const std::vector<FitResult>& results, const config_struct& cfg, TFile* o_file, const std::vector<std::string>& vals) {
    const auto& eta_bins = cfg.templ.eta_bins;
    const auto& pt_bins  = cfg.templ.pt_bins;

    if ((!o_file) || (o_file->IsZombie())) {
        std::cout << "ERROR: Invalid TFile pointer, exiting..." << std::endl;
        return 1;
    }

    std::string dir_name = "Plots_" + cfg.templ.bins_settup;
    
    if (cfg.analysis.sample_pass_data == "histo") {
        dir_name = dir_name + "_h";
    } else {
        dir_name = dir_name + "_d";
    }

    if (cfg.analysis.sample_pass_mc == "histo") {
        dir_name = dir_name + "_h";
    } else {
        dir_name = dir_name + "_d";
    }

    if (cfg.analysis.sample_fail_data == "histo") {
        dir_name = dir_name + "_h";
    } else {
        dir_name = dir_name + "_d";
    }

    if (cfg.analysis.sample_fail_mc == "histo") {
        dir_name = dir_name + "_h";
    } else {
        dir_name = dir_name + "_d";
    }
    
    TDirectory* bin_dir = o_file->GetDirectory(dir_name.c_str());
    
    if (!bin_dir) {
        bin_dir = o_file->mkdir(dir_name.c_str());
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
    
    o_file->Close();

    return 0;
}

// ------------------------------------------------------------------------------------------------------------------------------------

void SaveBinFitCanvas(RooRealVar& mll, RooCategory& sample, RooAbsData& data, RooSimultaneous& simPdf, const FitResult& res, TFile* o_file,
     const config_struct& cfg) {
    
    gROOT->SetBatch(kTRUE);
    
    if (!o_file || o_file->IsZombie()) {
        std::cerr << "ERROR: Invalid output TFile pointer in SaveBinFitCanvas!" << std::endl;
        return;
    }

    std::string dir_name = "Plots_" + cfg.templ.bins_settup;
    
    if (cfg.analysis.sample_pass_data == "histo") {
        dir_name = dir_name + "_h";
    } else {
        dir_name = dir_name + "_d";
    }

    if (cfg.analysis.sample_pass_mc == "histo") {
        dir_name = dir_name + "_h";
    } else {
        dir_name = dir_name + "_d";
    }

    if (cfg.analysis.sample_fail_data == "histo") {
        dir_name = dir_name + "_h";
    } else {
        dir_name = dir_name + "_d";
    }

    if (cfg.analysis.sample_fail_mc == "histo") {
        dir_name = dir_name + "_h";
    } else {
        dir_name = dir_name + "_d";
    }

    TDirectory* bin_dir = o_file->GetDirectory(dir_name.c_str());
    
    if (!bin_dir) {
        bin_dir = o_file->mkdir(dir_name.c_str());
    }

    bin_dir->cd();

    TCanvas c("c_fit", "Fit Pass and Fail", 1200, 600);
    c.Divide(2, 1);

    c.cd(1);
    gPad->SetLeftMargin(0.13);
   
    RooPlot* frame_pass = mll.frame(RooFit::Title("Category: PASS"));

    data.plotOn(frame_pass, RooFit::Cut("sample==sample::Pass"), RooFit::Name("data_pass"));
    simPdf.plotOn(frame_pass, RooFit::Slice(sample, "Pass"), RooFit::ProjWData(sample, data), RooFit::Name("pdf_pass"));
    simPdf.plotOn(frame_pass, RooFit::Slice(sample, "Pass"), 
              RooFit::Components("conv_pass_sig"), 
              RooFit::ProjWData(sample, data), 
              RooFit::LineColor(kGreen+2), RooFit::LineStyle(kDotted));
    simPdf.plotOn(frame_pass, RooFit::Slice(sample, "Pass"), RooFit::Components("bkg_pass"), RooFit::ProjWData(sample, data),
     RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
    
    double chi2_pass = frame_pass->chiSquare("pdf_pass", "data_pass", 5);
    frame_pass->Draw();

    TPaveText* txt_pass = new TPaveText(0.65, 0.65, 0.88, 0.88, "NDC");
    txt_pass->SetFillColor(0);
    txt_pass->SetTextAlign(12);
    txt_pass->AddText(Form("Status: %d", res.fit_status));
    txt_pass->AddText(Form("#chi^{2}/ndof: %.2f", chi2_pass));
    txt_pass->AddText(Form("Eff: %.3f #pm %.3f", res.efficiency, res.efficiency_err));
    txt_pass->Draw("same");

    c.cd(2);
    gPad->SetLeftMargin(0.13);
    RooPlot* frame_fail = mll.frame(RooFit::Title("Category: FAIL"));
    
    data.plotOn(frame_fail, RooFit::Cut("sample==sample::Fail"), RooFit::Name("data_fail"));
    simPdf.plotOn(frame_fail, RooFit::Slice(sample, "Fail"), RooFit::ProjWData(sample, data), RooFit::Name("pdf_fail"));
    simPdf.plotOn(frame_fail, RooFit::Slice(sample, "Fail"), 
              RooFit::Components("conv_fail_sig"), 
              RooFit::ProjWData(sample, data), 
              RooFit::LineColor(kGreen+2), RooFit::LineStyle(kDotted));
    simPdf.plotOn(frame_fail, RooFit::Slice(sample, "Fail"), RooFit::Components("bkg_fail"), RooFit::ProjWData(sample, data),
     RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));

    double chi2_fail = frame_fail->chiSquare("pdf_fail", "data_fail", 5);
    frame_fail->Draw();

    TPaveText* txt_fail = new TPaveText(0.65, 0.65, 0.88, 0.88, "NDC");    
    txt_fail->SetFillColor(0);
    txt_fail->SetTextAlign(12);
    txt_fail->AddText(Form("Status: %d", res.fit_status));
    txt_fail->AddText(Form("#chi^{2}/ndof: %.2f", chi2_fail));
    txt_fail->AddText(Form("#sigma: %.2f #pm %.2f GeV", res.sigma, res.sigma_err));
    txt_fail->Draw("same");

    std::string canvas_name = Form("fit_etaBin%d_ptBin%d", res.eta_bin_idx, res.pt_bin_idx);
    c.SetName(canvas_name.c_str());
    
    c.Write(canvas_name.c_str(), TObject::kOverwrite);

    delete frame_pass;
    delete frame_fail;
}