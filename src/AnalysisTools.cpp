#include "AnalysisTools.h"
#include "Utils.h"

#include "TH1D.h"

#include <map>
#include <memory>

#include <RooHistPdf.h>
#include <RooFormulaVar.h>

#include <RooGaussian.h>
#include <RooAddPdf.h>

#include <RooFFTConvPdf.h>
#include <RooKeysPdf.h>

#include <RooFitResult.h>
#include <RooMinimizer.h>

#include <RooPlot.h>
#include <RooHist.h>
#include <TPaveText.h>

#include <TCanvas.h>
#include <TLine.h>

// ------------------------------------------------------------------------------------------------------------------------------------

int BinnedTemplateMaker(ROOT::RDF::RNode node, const config_struct& cfg, const int dataset) {   
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

    // Names
    std::string name_pass = sample + "_h3_pass";
    std::string name_fail = sample + "_h3_fail";
    
    std::string name_entries_pass = sample + "_h2_entries_pass";
    std::string name_entries_fail = sample + "_h2_entries_fail";

    std::string title = "3D Histogram_" + sample;
    std::string title_entries = "2D Entries Histogram_" + sample;

    // Models
    ROOT::RDF::TH3DModel model_pass(name_pass.c_str(), title.c_str(), eta_bins.size() - 1, eta_bins.data(), pt_bins.size() - 1, pt_bins.data(),
     mll_bins.size() - 1, mll_bins.data() );
    
    ROOT::RDF::TH3DModel model_fail(name_fail.c_str(), title.c_str(), eta_bins.size() - 1, eta_bins.data(), pt_bins.size() - 1, pt_bins.data(),
     mll_bins.size() - 1, mll_bins.data() );

    ROOT::RDF::TH2DModel model_entries_pass(name_entries_pass.c_str(), title_entries.c_str(), eta_bins.size() - 1, eta_bins.data(), pt_bins.size() - 1, pt_bins.data());
    ROOT::RDF::TH2DModel model_entries_fail(name_entries_fail.c_str(), title_entries.c_str(), eta_bins.size() - 1, eta_bins.data(), pt_bins.size() - 1, pt_bins.data());
    
    // RNode
    node_hist = node_hist
        .Define(sample + "_Probe_Pt_Pass", sample + "_Probe_Pt[" + sample + "_Mask_Pass]")
        .Define(sample + "_Probe_Eta_Pass", sample + "_Probe_Eta[" + sample + "_Mask_Pass]")
        .Define(sample + "_Mll_Pass", sample + "_Mll[" + sample + "_Mask_Pass]")
        .Define(sample + "_Probe_Pt_Fail", sample + "_Probe_Pt[!" + sample + "_Mask_Pass]")
        .Define(sample + "_Probe_Eta_Fail", sample + "_Probe_Eta[!" + sample + "_Mask_Pass]")
        .Define(sample + "_Mll_Fail", sample + "_Mll[!" + sample + "_Mask_Pass]");

    if (cfg.templ.template_type.find("BINNED") != std::string::npos) {

        // Creating Smart Pointers
        auto h3_pass = node_hist.Histo3D(model_pass, sample + "_Probe_Eta_Pass", sample + "_Probe_Pt_Pass", sample + "_Mll_Pass");
        auto h3_fail = node_hist.Histo3D(model_fail, sample + "_Probe_Eta_Fail", sample + "_Probe_Pt_Fail", sample + "_Mll_Fail");
        
        auto h2_pass = node_hist.Histo2D(model_entries_pass, sample + "_Probe_Eta_Pass", sample + "_Probe_Pt_Pass");
        auto h2_fail = node_hist.Histo2D(model_entries_fail, sample + "_Probe_Eta_Fail", sample + "_Probe_Pt_Fail");

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

        o_template_file.Close();
    } else {
        std::cout << "Analysis DATA option selected, exiting from BinnedTemplateMaker." << std::endl;
        return 1; 
    }

    std::cout << "Saving binned data in the file " << cfg.templ.o_template_file_data << std::endl;

    return 0;
}

// ------------------------------------------------------------------------------------------------------------------------------------

int UnbinnedTemplateMaker(ROOT::RDF::RNode node, const config_struct& cfg) {
    ROOT::RDF::RNode node_roll = node;

    // Getting N cores
    int nMT = ROOT::GetThreadPoolSize(); 

    std::vector<std::vector<double>> pt(nMT);
    std::vector<std::vector<double>> eta(nMT);
    std::vector<std::vector<double>> mll(nMT);
    std::vector<std::vector<double>> pass_fail(nMT);
    
    for (int i = 0; i < nMT; ++i) {
        pt[i].reserve(100000);
        eta[i].reserve(100000);
        mll[i].reserve(100000);
        pass_fail[i].reserve(100000);
    }

    // Starts event loop
    node_roll.ForeachSlot([&pt, &eta, &mll, &pass_fail] (unsigned int slot, const ROOT::RVec<float>& RV_pt, const ROOT::RVec<float>& RV_eta,
            const ROOT::RVec<float>& RV_mll, const ROOT::RVec<bool>& RV_pass_fail) {
            
            for (size_t i = 0; i < RV_pt.size(); i++) {
                if (RV_pt[i] > 25) {
                    pt[slot].push_back(RV_pt[i]);
                    eta[slot].push_back(RV_eta[i]);
                    mll[slot].push_back(RV_mll[i]);
                    pass_fail[slot].push_back(RV_pass_fail[i]);
                }
            }
        }, {cfg.templ.dataset + "_Probe_Pt", cfg.templ.dataset + "_Probe_Eta", cfg.templ.dataset + "_Mll", cfg.templ.dataset + "_Mask_Pass"}
    );
    
    std::cout << "RDF phase done." << std::endl;

    size_t total_elements = 0;
    for (int i = 0; i < nMT; i++) {
        total_elements += pt[i].size();
    }
    
    std::cout << "Merging " << total_elements << " total particles from " << nMT << " threads..." << std::endl;

    std::vector<float> flat_pt;
    flat_pt.reserve(total_elements);
    
    std::vector<float> flat_eta;
    flat_eta.reserve(total_elements);
    
    std::vector<float> flat_mll;
    flat_mll.reserve(total_elements);
    
    std::vector<int> flat_pass_fail;
    flat_pass_fail.reserve(total_elements);

    // Loop in cores -> mergin data
    for (int i = 0; i < nMT; ++i) {
        flat_pt.insert(flat_pt.end(), std::make_move_iterator(pt[i].begin()), std::make_move_iterator(pt[i].end()));
        flat_eta.insert(flat_eta.end(), std::make_move_iterator(eta[i].begin()), std::make_move_iterator(eta[i].end()));
        flat_mll.insert(flat_mll.end(), std::make_move_iterator(mll[i].begin()), std::make_move_iterator(mll[i].end()));
        flat_pass_fail.insert(flat_pass_fail.end(), std::make_move_iterator(pass_fail[i].begin()), std::make_move_iterator(pass_fail[i].end()));
       
        // Security check
        pt[i].clear();
        pt[i].shrink_to_fit();

        eta[i].clear();
        eta[i].shrink_to_fit();

        mll[i].clear();
        mll[i].shrink_to_fit();

        pass_fail[i].clear();
        pass_fail[i].shrink_to_fit();
        
    }

    std::cout << "Merging done." << std::endl;

    TFile o_template_file(cfg.templ.o_template_file_data.c_str(), "UPDATE");

    if (o_template_file.IsZombie()) {
        std::cout << "ERROR: Cannot find output file: " << cfg.templ.o_template_file_data << std::endl;
        return 1;
    }

    o_template_file.cd();

    std::string tree_name = "Tree_" + cfg.templ.dataset + "_Flat";
    TTree tree(tree_name.c_str(), "Flattened tree");

    float b_pt, b_eta, b_mll;
    int b_pass_fail;

    tree.Branch("Probe_Pt", &b_pt, "Probe_Pt/F");
    tree.Branch("Probe_Eta", &b_eta, "Probe_Eta/F");
    tree.Branch("Probe_Mll", &b_mll, "Probe_Mll/F");
    tree.Branch("Mask_Pass", &b_pass_fail, "Mask_Pass/I");

    // Filling the flat tree
    for (size_t i = 0; i < total_elements; i++) {
        b_pt = flat_pt[i];
        b_eta = flat_eta[i];
        b_mll = flat_mll[i];
        b_pass_fail = flat_pass_fail[i];
        
        tree.Fill();
    }

    tree.Write("", TObject::kOverwrite);
    o_template_file.Close();// Flush
    
    std::cout << "Done. Saved " << total_elements << " entries." << std::endl;
    return 0;
}

// ------------------------------------------------------------------------------------------------------------------------------------
 
int LoadBinnedTemplate(const config_struct& cfg, std::vector<Template_RooF>& container) { 
    TH3D* h3_MC_pass{nullptr};
    TH3D* h3_MC_fail{nullptr};
    TH3D* h3_DATA_pass{nullptr};
    TH3D* h3_DATA_fail{nullptr};
 
    // Initializing LAZY flags
    bool h_mc_pass = false;
    bool h_mc_fail = false;
    bool h_data_pass = false;
    bool h_data_fail = false;
    
    bool d_mc_pass = false;
    bool d_mc_fail = false;
    bool d_data_pass = false;
    bool d_data_fail = false;
 
    if (cfg.analysis.sample_pass_mc == "binned") {
        h_mc_pass = true;
    } else if (cfg.analysis.sample_pass_mc == "unbinned") {
        d_mc_pass = true;
    }
    
    if (cfg.analysis.sample_pass_data == "binned") {
        h_data_pass = true;
    } else if (cfg.analysis.sample_pass_data == "unbinned") {
        d_data_pass = true;
    }
    
    if (cfg.analysis.sample_fail_mc == "binned") {
        h_mc_fail = true;
    } else if (cfg.analysis.sample_fail_mc == "unbinned") {
        d_mc_fail = true;
    }
    
    if (cfg.analysis.sample_fail_data == "binned") {
        h_data_fail = true;
    } else if (cfg.analysis.sample_fail_data == "unbinned") {
        d_data_fail = true;
    }
    
    // Opening file
    std::unique_ptr<TFile> file(TFile::Open(cfg.templ.o_template_file_data.c_str(), "READ"));
    
    if ((!file) || (file->IsZombie())) {
        std::cout << "ERROR: Cannot open the template input file: " << cfg.templ.o_template_file_data << std::endl;
        return 1;
    }
 
    TDirectory* bin_dir{nullptr};
    TTree* tree_data{nullptr};
    TTree* tree_mc{nullptr};
 
    // Directory of the binned data.
    if (h_data_pass || h_data_fail || h_mc_fail || h_mc_pass) {
        bin_dir = file->GetDirectory(cfg.templ.bins_settup.c_str());
        
        if (!bin_dir) {
            std::cout << "ERROR: read operation fails, exiting." << std::endl;
            return 1;
        }
    }

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
    
    // Unbinned data
    if (d_data_pass || d_data_fail) {
        std::string tree_name = "Tree_DATA_Flat";
        tree_data = file->Get<TTree>(tree_name.c_str());
    }

    if (d_mc_pass || d_mc_fail) {
        std::string tree_name = "Tree_MC_Flat";
        tree_mc = file->Get<TTree>(tree_name.c_str());
    }

    // Starting bin settup
    std::vector<float> pt_bins = cfg.templ.pt_bins;
    std::vector<float> eta_bins = cfg.templ.eta_bins;
    
    const int n_pt_bins = pt_bins.size() - 1;
    const int n_eta_bins = eta_bins.size() - 1;
    const int tot_bins = n_pt_bins * n_eta_bins;
 
    container.clear();
    container.reserve(tot_bins);
 
    std::cout << "Starting load ..." << std::endl;
 
    for (int i = 0; i < n_pt_bins; i++) {
        int bin_pt = i + 1;
 
        for (int j = 0; j < n_eta_bins; j++) {
            int bin_eta = j + 1;
 
            // For memory control
            std::unique_ptr<TH1D> h1_mc_pass;
            std::unique_ptr<TH1D> h1_mc_fail;
            std::unique_ptr<TH1D> h1_data_pass;
            std::unique_ptr<TH1D> h1_data_fail;
 
            if (h_mc_pass) {
                std::string name = "h_mll_mc_pass_eta" + std::to_string(bin_eta) + "_pt" + std::to_string(bin_pt);
                TH1D* hist = h3_MC_pass->ProjectionZ(name.c_str(), bin_eta, bin_eta, bin_pt, bin_pt);
                
                hist->SetDirectory(nullptr);
                
                h1_mc_pass.reset(hist);
            }

            if (h_mc_fail) {
                std::string name = "h_mll_mc_fail_eta" + std::to_string(bin_eta) + "_pt" + std::to_string(bin_pt);
                TH1D* hist = h3_MC_fail->ProjectionZ(name.c_str(), bin_eta, bin_eta, bin_pt, bin_pt);
                
                hist->SetDirectory(nullptr);
                
                h1_mc_fail.reset(hist);
            }
            
            if (h_data_pass) {
                std::string name = "h_mll_data_pass_eta" + std::to_string(bin_eta) + "_pt" + std::to_string(bin_pt);
                TH1D* hist = h3_DATA_pass->ProjectionZ(name.c_str(), bin_eta, bin_eta, bin_pt, bin_pt);
                
                hist->SetDirectory(nullptr);
                
                h1_data_pass.reset(hist);
            }
            
            if (h_data_fail) {
                std::string name = "h_mll_data_fail_eta" + std::to_string(bin_eta) + "_pt" + std::to_string(bin_pt);
                TH1D* hist = h3_DATA_fail->ProjectionZ(name.c_str(), bin_eta, bin_eta, bin_pt, bin_pt);
                
                hist->SetDirectory(nullptr);
                
                h1_data_fail.reset(hist);
            }
 
            container.emplace_back(Template_RooF{
                bin_eta,
                bin_pt,
                std::move(h1_mc_pass),
                std::move(h1_mc_fail),
                std::move(h1_data_pass),
                std::move(h1_data_fail),
                nullptr,
                nullptr,
                nullptr,
                nullptr
            });
        }
    }
  
    std::cout << "Loading unbinned data..." << std::endl;

    if (d_mc_pass || d_mc_fail) {
        LoadUnbinnedTemplate(tree_mc, cfg, 2, container);
    }
    if (d_data_pass || d_data_fail) {
        LoadUnbinnedTemplate(tree_data, cfg, 1, container);
    }

    std::cout << "Load finished, exiting from LoadBinnedTemplate..." << std::endl;
 
    file->Close();
    return 0;
}

// ------------------------------------------------------------------------------------------------------------------------------------

int LoadUnbinnedTemplate(TTree* tree, const config_struct& cfg, const int dataset, std::vector<Template_RooF>& container) {
    if (!tree) {
        std::cout << "ERROR: tree invalid pointer, exiting..." << std::endl;
        return 1;
    }

    const std::vector<float> pt_bins = cfg.templ.pt_bins;
    const std::vector<float> eta_bins = cfg.templ.eta_bins;

    float min_pt = pt_bins.front();
    float max_pt = pt_bins.back();

    float min_eta = eta_bins.front();
    float max_eta = eta_bins.back();

    bool process_pass = false;
    bool process_fail = false;

    if (dataset == 1) {
        process_pass = (cfg.analysis.sample_pass_data == "unbinned");
        process_fail = (cfg.analysis.sample_fail_data == "unbinned");
    }

    if (dataset == 2) {
        process_pass = (cfg.analysis.sample_pass_mc == "unbinned");
        process_fail = (cfg.analysis.sample_fail_mc == "unbinned");
    }  

    int n_bins = (pt_bins.size() - 1) * (eta_bins.size() - 1);
    int n_final_bins = 2 * n_bins;

    std::vector<std::vector<float>> mll_buffers(n_final_bins);
    
    if (process_fail) {
        for (int i = 0; i < n_bins; i++) {
            mll_buffers[i].reserve(20000);
        }
    }

    if (process_pass) {
        for (int i = n_bins; i < n_final_bins; i++) {
            mll_buffers[i].reserve(20000);
        }
    }

    float t_pt, t_eta, t_mll;
    int t_mask;

    tree->SetBranchAddress("Probe_Pt", &t_pt);
    tree->SetBranchAddress("Probe_Eta", &t_eta);
    tree->SetBranchAddress("Probe_Mll", &t_mll);
    tree->SetBranchAddress("Mask_Pass", &t_mask);

    Long64_t nentries = tree->GetEntries();
    for (Long64_t entry = 0; entry < nentries; entry++) {
        tree->GetEntry(entry);

        if (t_pt < min_pt || t_pt >= max_pt || t_eta < min_eta || t_eta >= max_eta) {
            continue;
        }

        bool is_pass = (t_mask == 1);

        if ((is_pass) && (!process_pass)) {
            continue;
        }
        if ((!is_pass) && (!process_fail)) {
            continue;
        }
        
        int i_pt = std::distance(pt_bins.begin(),  std::upper_bound(pt_bins.begin(), pt_bins.end(), t_pt)) - 1;
        int i_eta = std::distance(eta_bins.begin(), std::upper_bound(eta_bins.begin(), eta_bins.end(), t_eta)) - 1;

        int k = i_eta + (i_pt * (eta_bins.size() - 1)) + (t_mask * (pt_bins.size() - 1) * (eta_bins.size() - 1));

        mll_buffers[k].push_back(t_mll);
    }

    RooRealVar mll("mll", "Invariant Mass", 60.0, 120.0);
    RooArgSet vars(mll);

    for (int i = 0; i < n_final_bins; i++) {
        if (mll_buffers[i].empty()) {
            continue;
        }

        int idx = i % n_bins;
        bool is_pass = (i >= n_bins);

        std::string sample_tag = (dataset == 1) ? "DATA" : "MC";
        std::string pass_tag = (is_pass) ? "Pass" : "Fail";
        std::string data_name = "d_" + sample_tag + "_" + pass_tag + "_idx" + std::to_string(idx);

        auto data = std::make_unique<RooDataSet>(data_name.c_str(), data_name.c_str(), vars);

        for (float mass : mll_buffers[i]) {
            mll.setVal(mass);
            data->add(vars);
        }

        // DATA
        if (dataset == 1) { 
            if (is_pass) {
                container[idx].d_DATA_pass = std::move(data);
            } else {
                container[idx].d_DATA_fail = std::move(data);
            }
        // MC
        } else if (dataset == 2) {
            if (is_pass) {
                container[idx].d_MC_pass = std::move(data);
            } else {
                container[idx].d_MC_fail = std::move(data);
            }
        }
    }
    return 0;
}

// ------------------------------------------------------------------------------------------------------------------------------------

int SaveMapFittedValues(TFile* o_file, const std::vector<FitResult>& results, const config_struct& cfg, const std::vector<std::string>& vals) {
    const auto& eta_bins = cfg.templ.eta_bins;
    const auto& pt_bins  = cfg.templ.pt_bins;

    if ((!o_file) || (o_file->IsZombie())) {
        std::cout << "ERROR: Invalid TFile pointer, exiting..." << std::endl;
        return 1;
    }

    std::string dir_name = "Plots_" + cfg.templ.bins_settup;
    
    if (cfg.analysis.sample_pass_data == "binned") {
        dir_name = dir_name + "_h";
    } else {
        dir_name = dir_name + "_d";
    }

    if (cfg.analysis.sample_pass_mc == "binned") {
        dir_name = dir_name + "_h";
    } else {
        dir_name = dir_name + "_d";
    }

    if (cfg.analysis.sample_fail_data == "binned") {
        dir_name = dir_name + "_h";
    } else {
        dir_name = dir_name + "_d";
    }

    if (cfg.analysis.sample_fail_mc == "binned") {
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

            if (var_names == "fit_status") {
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

    if (!o_file || o_file->IsZombie()) {
        std::cout << "ERROR: Invalid output TFile pointer in SaveBinFitCanvas!" << std::endl;
        return;
    }

    std::string dir_name = "Plots_" + cfg.templ.bins_settup;
    
    if (cfg.analysis.sample_pass_data == "binned") {
        dir_name = dir_name + "_h";
    } else {
        dir_name = dir_name + "_d";
    }

    if (cfg.analysis.sample_pass_mc == "binned") {
        dir_name = dir_name + "_h";
    } else {
        dir_name = dir_name + "_d";
    }

    if (cfg.analysis.sample_fail_data == "binned") {
        dir_name = dir_name + "_h";
    } else {
        dir_name = dir_name + "_d";
    }

    if (cfg.analysis.sample_fail_mc == "binned") {
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

    // -----
    // Pad 1
    // -----

    c.cd(1);
    gPad->SetLeftMargin(0.13);
   
    std::unique_ptr<RooPlot> frame_pass(mll.frame(RooFit::Title("Category: PASS")));

    data.plotOn(
        frame_pass.get(), 
        RooFit::Cut("sample==sample::Pass"), 
        RooFit::Name("data_pass"));
    
    simPdf.plotOn(
        frame_pass.get(), 
        RooFit::Slice(sample, "Pass"), 
        RooFit::ProjWData(sample, data), 
        RooFit::Name("pdf_pass"));
    
    simPdf.plotOn(
        frame_pass.get(), 
        RooFit::Slice(sample, "Pass"), 
        RooFit::Components("conv_pass_sig"), 
        RooFit::ProjWData(sample, data), 
        RooFit::LineColor(kGreen + 2), RooFit::LineStyle(kDotted));
    
    simPdf.plotOn(
        frame_pass.get(), 
        RooFit::Slice(sample, "Pass"), 
        RooFit::Components("bkg_pass"), 
        RooFit::ProjWData(sample, data),
        RooFit::LineStyle(kDashed), 
        RooFit::LineColor(kRed));
    
    
    double chi2_pass = frame_pass->chiSquare("pdf_pass", "data_pass", 5);
    
    frame_pass->Draw();

    TPaveText* txt_pass = new TPaveText(0.65, 0.65, 0.88, 0.88, "NDC");
    txt_pass->SetFillColor(0);
    txt_pass->SetTextAlign(12);
    txt_pass->AddText(Form("Status: %d", res.fit_status));
    txt_pass->AddText(Form("#chi^{2}/ndof: %.2f", chi2_pass));
    txt_pass->AddText(Form("Eff: %.3f #pm %.3f", res.efficiency, res.efficiency_err));
    txt_pass->Draw("same");

    // -----
    // Pad 2
    // -----

    c.cd(2);
    gPad->SetLeftMargin(0.13);
    
    std::unique_ptr<RooPlot> frame_fail(mll.frame(RooFit::Title("Category: FAIL")));
    
    data.plotOn(
        frame_fail.get(), 
        RooFit::Cut("sample==sample::Fail"), 
        RooFit::Name("data_fail"));
    
    simPdf.plotOn(
        frame_fail.get(),
        RooFit::Slice(sample, "Fail"), 
        RooFit::ProjWData(sample, data), 
        RooFit::Name("pdf_fail"));
    
    simPdf.plotOn(
        frame_fail.get(), 
        RooFit::Slice(sample, "Fail"), 
        RooFit::Components("conv_fail_sig"), 
        RooFit::ProjWData(sample, data), 
        RooFit::LineColor(kGreen + 2), 
        RooFit::LineStyle(kDotted));
    
    simPdf.plotOn(
        frame_fail.get(),
        RooFit::Slice(sample, "Fail"), 
        RooFit::Components("bkg_fail"), 
        RooFit::ProjWData(sample, data),
        RooFit::LineStyle(kDashed), 
        RooFit::LineColor(kRed));

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
}

// ------------------------------------------------------------------------------------------------------------------------------------

int EfficiencyFitter(std::vector<Template_RooF>& analysis_struct, const config_struct& cfg, std::vector<FitResult>& results, TFile* o_file) {
    
    if (!o_file) {
        std::cout << "ERROR: invalid output file, exiting..." << std::endl;
        return 0;
    }

    // -----
    // Start
    // -----

    results.reserve(analysis_struct.size());
    
    RooRealVar mll("mll", "m_{#mu+#mu-}", 60.0, 120.0, "GeV");
    mll.setBins(8192, "cache");// For FFT calculus

    RooCategory sample("sample", "TagAndProbe categories");
    sample.defineType("Pass");
    sample.defineType("Fail");

    int successful_fits = 0;
    int failed_fits = 0;

    // ------------------------------
    // Loop on P_t_Z0 and Eta_Z0 bins
    // ------------------------------

    for (const auto& it : analysis_struct) {
        std::cout << "Fitting -> Eta bin: " << it.eta_bin_idx << " , Pt bin: " << it.pt_bin_idx << std::endl;

        // ------------------------------------
        // MC PASS option -> BINNED or UNBINNED
        // ------------------------------------

        std::unique_ptr<RooAbsPdf> mc_pass_pdf{nullptr};
        std::unique_ptr<RooDataHist> hist_mc_pass{nullptr};

        if (cfg.analysis.sample_pass_mc == "unbinned") {
            mc_pass_pdf = std::make_unique<RooKeysPdf>("mc_pass_pdf", "MC Pass KeysPdf", mll, *it.d_MC_pass, RooKeysPdf::NoMirror, 1.5);

            if (!mc_pass_pdf) {
                std::cout << "ERROR: invalid mc_pass_pdf pointer, exiting..." << std::endl; 
                return 1;
            }

        } else if (cfg.analysis.sample_pass_mc == "binned") {
            hist_mc_pass = std::make_unique<RooDataHist>("hist_mc_pass", "Pass data histogram", mll, it.h_MC_pass.get());
            auto hist_pdf = std::make_unique<RooHistPdf>("hist_mc_pass_pdf", "Pass data pdf", mll, *hist_mc_pass);
            
            hist_pdf->setInterpolationOrder(2);
            mc_pass_pdf = std::move(hist_pdf);
            
            if (!mc_pass_pdf) {
                std::cout << "ERROR: invalid mc_pass_pdf pointer, exiting..." << std::endl; 
                return 1;
            }
        }

        // ------------------------------------
        // MC FAIL option -> BINNED or UNBINNED
        // ------------------------------------

        std::unique_ptr<RooAbsPdf> mc_fail_pdf = nullptr;
        std::unique_ptr<RooDataHist> hist_mc_fail = nullptr;

        if (cfg.analysis.sample_fail_mc == "unbinned") {            
            mc_fail_pdf = std::make_unique<RooKeysPdf>("mc_fail_pdf", "MC Fail KeysPdf", mll, *it.d_MC_fail, RooKeysPdf::NoMirror, 1.5);

            if (!mc_fail_pdf) {
                std::cout << "ERROR: invalid mc_fail_pdf pointer, exiting..." << std::endl; 
                return 1;
            }
        } else if (cfg.analysis.sample_fail_mc == "binned") {
            hist_mc_fail = std::make_unique<RooDataHist>("hist_mc_fail", "Pass fail histogram", mll, it.h_MC_fail.get());
            auto hist_pdf = std::make_unique<RooHistPdf>("hist_mc_fail_pdf", "Pass fail pdf", mll, *hist_mc_fail);
            
            hist_pdf->setInterpolationOrder(2);
            mc_fail_pdf = std::move(hist_pdf);

            if (!mc_fail_pdf) {
                std::cout << "ERROR: invalid mc_fail_pdf pointer, exiting..." << std::endl; 
                return 1;
            }
        }

        // -------------------------------------------
        // DATA -> PASSED/FAILED -> SAME SAMPLE OPTION
        // -------------------------------------------

        std::unique_ptr<RooAbsData> sig_data;

        if ((cfg.analysis.sample_pass_data == "unbinned") && (cfg.analysis.sample_fail_data == "unbinned")) {    
            sig_data = std::make_unique<RooDataSet>("sig_data", "Signal data unbinned", 
                RooArgSet(mll),
                RooFit::Index(sample),
                RooFit::Import("Pass", *it.d_DATA_pass), 
                RooFit::Import("Fail", *it.d_DATA_fail));
            
            if (!sig_data) {
                std::cout << "ERROR: invalid sig_data pointer, exiting..." << std::endl; 
                return 1;
            }
            
        } else if ((cfg.analysis.sample_pass_data == "binned") && (cfg.analysis.sample_fail_data == "binned")){
            std::cout << "Data binned loaded" << std::endl;
        
            std::map<std::string, TH1*> map_hist_data;
            map_hist_data["Pass"] = it.h_DATA_pass.get();
            map_hist_data["Fail"] = it.h_DATA_fail.get();

            sig_data = std::make_unique<RooDataHist>("sig_data", "Signal data binned", mll, sample, map_hist_data);
            
            if (!sig_data) {
                std::cout << "ERROR: invalid sig_data pointer, exiting..." << std::endl; 
                return 1;
            }
        } else {
            std::cout << "ERROR: invalid combination in JSON file, exiting..." << std::endl;
            return 1;
        }

        // -----------------------------------------
        // Getting INITIAL value for YIELD parameter
        // -----------------------------------------
         
        bool is_data_unbinned = (cfg.analysis.sample_pass_data == "unbinned");

        double data_pass_entries = (is_data_unbinned) ? it.d_DATA_pass->sumEntries() : it.h_DATA_pass->Integral();
        double data_fail_entries = (is_data_unbinned) ? it.d_DATA_fail->sumEntries() : it.h_DATA_fail->Integral();
        double total_entries = data_pass_entries + data_fail_entries;

        // -----------------------------
        // DEFINING FIT MODEL/PARAMETERS
        // -----------------------------

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

        // --------------
        // PRE-FIT OPTION 
        // --------------

        if (cfg.analysis.pre_fit) {
            std::cout << "Starting pre-fitting." << std::endl;

            // Fixing efficiency value
            efficiency.setVal(1.0);
            efficiency.setConstant(kTRUE);
            RooDataHist hist_data_pass("hist_data_pass", "Pass data histogram", mll, it.h_DATA_pass.get());

            std::unique_ptr<RooAbsReal> nll_pass(model_pass.createNLL(
                hist_data_pass,
                RooFit::Extended(kTRUE), 
                RooFit::NumCPU(8)
            ));

            RooMinimizer m_pass(*nll_pass);

            m_pass.setStrategy(2);
            m_pass.setEps(1e-2);
            int status_migrad = m_pass.migrad();
        
            if (status_migrad == 0 || status_migrad == 1) {
                m_pass.hesse();
            }

            std::unique_ptr<RooFitResult> prefitRes(m_pass.save());

            // Passing fitted values as parameter inputs -> mean/sigma convolutional gaussian
            double mu_fit = mu.getVal();
            double mu_err = mu.getError();

            double sigma_fit = sigma.getVal();
            double sigma_err = sigma.getError();

            mu.setConstant(kFALSE);
            sigma.setConstant(kFALSE);

            mu.setVal(mu_fit);
            sigma.setVal(sigma_fit);
        
            // Unfixing efficiency
            efficiency.setConstant(kFALSE);
            efficiency.setVal(cfg.analysis.params.efficiency.at(0));
        }

        // -----------------------
        // FINAL FIT -> NLL OPTION
        // -----------------------

        std::unique_ptr<RooAbsReal> nll(simPdf.createNLL(
            *sig_data,
            RooFit::Extended(kTRUE), 
            RooFit::NumCPU(8)
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

        // Verbose
        raw_fitRes->Print("v");

        // ---------------------
        // FILLING RESULT STRUCT
        // ---------------------

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

        // ---------------------------------
        // SAVING FIT PLOTS -> QUALITY CHECK
        // ---------------------------------

        SaveBinFitCanvas(mll, sample, *sig_data, simPdf, results.back(), o_file, cfg);
    }
    
    std::cout << "Succeded fits: " << successful_fits << std::endl;
    std::cout << "Failed fits: " << failed_fits << std::endl;

    std::cout << "Fit procedure finished." << std::endl;
    
    return 0;
}

// ------------------------------------------------------------------------------------------------------------------------------------

void SaveEventFitCanvas(RooRealVar& mll, RooAbsPdf& model, RooAbsData& data, RooAbsPdf& bkg_pdf, const EventFitResult& res, TH1D* h_mll,
TDirectory* o_dir, const int bin_idx, const std::string& tag) {
    
    if (!o_dir) {
        return;
    }

    std::cout << "Saving histogram for bin: " << bin_idx << std::endl;
    o_dir->cd();

    int n_bins_frame = h_mll->GetNbinsX();

    TCanvas canvas(("canvas_" + tag + "_bin_" + std::to_string(bin_idx)).c_str(), "", 800, 800);
    
    TPad* pad_main  = new TPad("pad_main", "", 0.0, 0.30, 1.0, 1.0);
    TPad* pad_resid = new TPad("pad_resid", "", 0.0, 0.0, 1.0, 0.30);
    
    pad_main->SetBottomMargin(0.02);
    pad_resid->SetTopMargin(0.02);
    pad_resid->SetBottomMargin(0.35);
    pad_main->Draw();
    pad_resid->Draw();

    // ----------------
    // MAIN PAD (Fit)
    // ----------------
    pad_main->cd();

    RooPlot* frame = mll.frame(
        RooFit::Title(("mll_" + tag + "_bin_" + std::to_string(bin_idx)).c_str()), 
        RooFit::Bins(n_bins_frame)
    );

    data.plotOn(frame, RooFit::Name("data_hist"));
    model.plotOn(frame, RooFit::Name("model"), RooFit::LineColor(kBlue));
    model.plotOn(frame, RooFit::Components(bkg_pdf), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));

    frame->GetXaxis()->SetLabelSize(0);
    frame->GetXaxis()->SetTitle("");
    frame->Draw();

    // ---------------------
    // RESIDUALS PAD (Pulls)
    // ---------------------
    pad_resid->cd();

    RooHist* hpull = frame->pullHist("data_hist", "model");
    
    RooPlot* frame_pull = mll.frame(RooFit::Title(""));
    
    if (hpull) {
        frame_pull->addPlotable(hpull, "P");
    }

    frame_pull->GetYaxis()->SetTitle("Pull");
    frame_pull->GetYaxis()->SetTitleSize(0.11);
    frame_pull->GetYaxis()->SetLabelSize(0.09);
    frame_pull->GetYaxis()->SetTitleOffset(0.45);
    frame_pull->GetYaxis()->SetNdivisions(504);
    frame_pull->GetXaxis()->SetTitleSize(0.12);
    frame_pull->GetXaxis()->SetLabelSize(0.10);
    frame_pull->GetYaxis()->SetRangeUser(-5.0, 5.0);
    frame_pull->Draw();

    TLine* zero_line = new TLine(mll.getMin(), 0.0, mll.getMax(), 0.0);
    zero_line->SetLineColor(kGray+2);
    zero_line->SetLineStyle(2);
    zero_line->Draw("same");

    canvas.Write(("canvas_" + tag + "_bin_" + std::to_string(bin_idx)).c_str(), TObject::kOverwrite);

    delete frame;
    delete frame_pull;
}

// ------------------------------------------------------------------------------------------------------------------------------------

std::vector<std::unique_ptr<TH1D>> PrepareEventFit(EventHisto& event_histo, const std::string& tag) {

    // Binned data container -> input
    std::vector<std::unique_ptr<TH1D>> container;
    
    // Central dataset [ X(P_t_Z0/Y_Z0/Phi*_Z0) , Y(Invariant Mass) ]
    TH2D* histo;

    // Kinematical quantity selection.
    if (tag == "pt") {
        histo = event_histo.h2_mll_pt.GetPtr();
    } else if (tag == "y") {
        histo = event_histo.h2_mll_y.GetPtr();
    } else if (tag == "phis") {
        histo = event_histo.h2_mll_phis.GetPtr();
    } else {
        std::cout << "ERROR: invalid tag input, exiting..." << std::endl;
    }

    int n_bins = histo->GetNbinsX();
    
    container.reserve(n_bins);

    // Saving Proyections of Invariant Mass (mll)
    for (int i = 0; i < n_bins; i++) {
        int bin_root = i + 1;

        std::string name = "h_mll_" + tag + "_bin" + std::to_string(i);
        TH1D* histo_mll = histo->ProjectionY(name.c_str(), bin_root, bin_root);

        histo_mll->SetDirectory(nullptr);

        container.emplace_back(std::unique_ptr<TH1D>(histo_mll));
    }

    return container;
}

// ------------------------------------------------------------------------------------------------------------------------------------

std::vector<std::unique_ptr<RooDataSet>> PrepareEventFitModel(TTree* tree, const config_struct& cfg, const std::string& tag) {
    std::vector<std::unique_ptr<RooDataSet>> container;
    std::vector<double> vector_bins;

    if (cfg.unfold.use_custom_bins == true) {
        if (tag == "pt") {
            vector_bins = cfg.unfold.pt_bins.reco_vec;
        } else if (tag == "y") {
            vector_bins = cfg.unfold.y_bins.reco_vec;
        } else if (tag == "phis") {
            vector_bins = cfg.unfold.phis_bins.reco_vec;
        }
    } else {
        if (tag == "pt") {
            const auto& pt = cfg.unfold.pt_bins;
            vector_bins = CreateBins(pt.reco_bins, pt.min, pt.max, pt.distribution, pt.split);

        } else if (tag == "y") {
            const auto& y = cfg.unfold.y_bins;
            vector_bins = CreateBins(y.reco_bins, y.min, y.max, y.distribution, y.split);            

        } else if (tag == "phis") {
            const auto& phis = cfg.unfold.phis_bins;
            vector_bins = CreateBins(phis.reco_bins, phis.min, phis.max, phis.distribution, phis.split);
        }
    }

    float min_pt = vector_bins.front();
    float max_pt = vector_bins.back();

    int n_bins = vector_bins.size() - 1;
    
    container.reserve(n_bins);

    std::vector<std::vector<double>> mll_buffers(n_bins);
    
    
    for (int i = 0; i < n_bins; i++) {
        mll_buffers[i].reserve(20000);
    }

    float t_val, t_mll;

    tree->SetBranchAddress("InvariantMass", &t_mll);
    
    if (tag == "pt") {
        tree->SetBranchAddress("Pt_Z", &t_val);
    } else if (tag == "y") {
        tree->SetBranchAddress("Y_Z", &t_val);
    } else if (tag == "phis") {
        tree->SetBranchAddress("Phis_Z", &t_val);
    }

    Long64_t nentries = tree->GetEntries();
    for (Long64_t entry = 0; entry < nentries; entry++) {
        tree->GetEntry(entry);
        
        int i_vector = std::distance(vector_bins.begin(),  std::upper_bound(vector_bins.begin(), vector_bins.end(), t_val)) - 1;
        if ((i_vector < 0) || (i_vector >= n_bins)) {
           continue;
        }

        mll_buffers[i_vector].push_back(t_mll);
    }

    RooRealVar mll("mll", "Invariant Mass", 60.0, 120.0);
    RooArgSet vars(mll);

    for (int i = 0; i < n_bins; i++) {
        if (mll_buffers[i].empty()) {
            continue;
        }

        std::string sample_tag = "MC";
        std::string data_name = "d_" + sample_tag + "_idx" + std::to_string(i);

        auto data = std::make_unique<RooDataSet>(data_name.c_str(), data_name.c_str(), vars);

        for (float mass : mll_buffers[i]) {
            mll.setVal(mass);
            data->add(vars);
        }

        container.push_back(std::move(data));
    }

    return container;
}

// ------------------------------------------------------------------------------------------------------------------------------------

EventFitResult EventSingleFit(int bin_idx, TH1D* h_mll, RooDataSet* d_mll_model, TDirectory* o_dir, const std::string& tag, const bool save_plots) {
    EventFitResult result;

    result.bin_idx = bin_idx;
    
    // Unpacking histogram
    if (!h_mll || !d_mll_model) {
        std::cout << "ERROR: invalid input histogram, exiting..." << std::endl;
        return result;
    }

    // Defining variable
    RooRealVar mll("mll", "m_{#mu#mu}", 60.0, 120.0, "GeV");

    // Defining data
    RooDataHist data_hist(("h_mll_" + tag + "_bin_" + std::to_string(bin_idx)).c_str(), "", mll, h_mll);

    // Signal distribution
    RooKeysPdf sig_pdf("sig_pdf", "Signal PDF from MC (KDE)", mll, *d_mll_model, RooKeysPdf::MirrorBoth);
    
    //RooHistPdf sig_pdf("sig_pdf", "Signal PDF from MC", mll, hist_mc_model);
    //sig_pdf.setInterpolationOrder(3);    
    
    // Bkg distribution -> exp
    RooRealVar lambda("lambda", "Decay constant", -0.02, -0.5, 0.0);
    RooExponential bkg_pdf("bkg_pdf", "Exponential Bkg", mll, lambda);

    double total_entries = h_mll->Integral();
    // Yields
    RooRealVar n_sig("n_sig", "Signal Yield", total_entries * 0.9, 0.0, total_entries * 1.5);
    RooRealVar n_bkg("n_bkg", "Bkg Yield", total_entries * 0.1, 0.0, total_entries * 1.5);

    // MODEL
    RooAddPdf model("model", "Signal + Bkg", RooArgList(sig_pdf, bkg_pdf), RooArgList(n_sig, n_bkg));

    // Fit
    std::unique_ptr<RooFitResult> fit_res(model.fitTo(
        data_hist, 
        RooFit::Extended(true),
        RooFit::NumCPU(8), 
        RooFit::Save(true)
    ));

    // Saving values
    result.n_sig = n_sig.getVal();
    result.n_sig_err = n_sig.getError();
    
    result.n_bkg = n_bkg.getVal();
    result.n_bkg_err = n_bkg.getError();
    
    result.lambda = lambda.getVal();
    result.lambda_err = lambda.getError();
    
    result.fit_status = fit_res ? fit_res->status() : -1;

    if (save_plots) {
        SaveEventFitCanvas(mll, model, data_hist, bkg_pdf, result, h_mll, o_dir, bin_idx, tag);
    }
    
    return result;
}

// ------------------------------------------------------------------------------------------------------------------------------------

std::vector<EventFitResult> EventFitWrapper(std::vector<std::unique_ptr<TH1D>>& container, std::vector<std::unique_ptr<RooDataSet>>& container_model, const std::string& tag, const bool save_plots, TDirectory* o_dir) {
    // Fit global container
    std::vector<EventFitResult> results;

    results.reserve(container.size());

    int succed_fits = 0;
    int failed_fits = 0;

    // Loop on Input histograms -> check needed(?)
    for (int i = 0; i < container.size(); i++) {
        
        // Getting smart pointer
        TH1D* histo = container[i].get();
        RooDataSet* data_model = container_model[i].get();

        // Fit
        EventFitResult res = EventSingleFit(i, histo, data_model, o_dir, tag, save_plots);

        // Fit status + counting
        if (res.fit_status == 0) {
            succed_fits++;
        } else {
            failed_fits++;
        }

        results.push_back(res);
    }

    std::cout << "Fits succeded: " << succed_fits << ", failed: " << failed_fits << std::endl;
    
    return results;
}