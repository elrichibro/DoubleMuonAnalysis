#include "Manager.h"

#include <string>
#include <vector>

#include <chrono>

// ------------------------------------------------------------------------------------------------------------------------------------
// ObjectTH1 Methods
// ------------------------------------------------------------------------------------------------------------------------------------

void ObjectTH1::Write(TFile& file) {
    if (th1) {
        th1->Write();
    }
}

void ObjectTH1::Draw(TCanvas& canvas) {
    if (!th1) {
        return;
    }
    canvas.cd();
    th1->Draw("E HIST");
}

// ------------------------------------------------------------------------------------------------------------------------------------
// ObjectTH2 Methods
// ------------------------------------------------------------------------------------------------------------------------------------

void ObjectTH2::Write(TFile& file) {
    if (th2) {
        th2->Write();
    }
}

void ObjectTH2::Draw(TCanvas& canvas) {
    if (!th2) {
        return;
    }
    canvas.cd();
    canvas.SetRightMargin(0.14); 
    th2->Draw("COLZ");
}

// ------------------------------------------------------------------------------------------------------------------------------------
// OutputManager Methods
// ------------------------------------------------------------------------------------------------------------------------------------

void OutputManager::AddToPipeline(const std::string& name, ROOT::RDF::RResultPtr<TH1D> hist) {
    pipeline.push_back(std::make_unique<ObjectTH1>(name, hist));
}

void OutputManager::AddToPipeline(const std::string& name, ROOT::RDF::RResultPtr<TH2D> hist) {
    pipeline.push_back(std::make_unique<ObjectTH2>(name, hist));
}

void OutputManager::Run() {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    for (auto& snap : snapshot_vec) {
        snap.GetValue(); 
        
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> DeltaT = end_time - start_time;
        
        std::cout << "Snapshot time: " << DeltaT.count() << std::endl;
    }

    std::unique_ptr<TFile> file_plots = nullptr;
    if (save_sel_plots) {
        file_plots = std::make_unique<TFile>(o_file_plots.c_str(), "RECREATE");
        std::cout << "Saving plots in file: " <<  o_file_plots << std::endl;
    }

    int i = 0;

    for (auto& it : pipeline) {
        it->Process(); 

        // Saving
        if (save_sel_plots) {
            it->Write(*file_plots); 
        }

        // Visualization
        if (visualize) {
            std::string c_name = "c_" + it->GetName();
            TCanvas* vis_canvas = new TCanvas(c_name.c_str(), it->GetName().c_str(), 800, 600);
            
            it->Draw(*vis_canvas);
            vis_canvas->Update();
        }
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> DeltaT = end_time - start_time;

        std::cout << "PipelineObj number: " << i << ", execution time: " << DeltaT.count() << std::endl;
        i++;
    }

    if (file_plots && file_plots->IsOpen()) {
        file_plots->Close();
        
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> DeltaT = end_time - start_time;
        
        std::cout << "Close file time: " << DeltaT.count() << std::endl;
    }
    
    std::cout << "Ending of Run" << std::endl;
}

void OutputManager::BookAnalysis(ROOT::RDF::RNode node, const config_struct& cfg) {
    
    std::string title_pt = ";" + cfg.pt_plot.title_axis + ";Efficiency;";
    std::string name_pt = cfg.general.dataset + "_p_{T}";
    ROOT::RDF::TH1DModel model_1D_pt(name_pt.c_str(), title_pt.c_str(), cfg.pt_plot.nbins, cfg.pt_plot.axis_min, 
    cfg.pt_plot.axis_max);

    std::string title_eta = ";" + cfg.eta_plot.title_axis + ";Efficiency;";
    std::string name_eta = cfg.general.dataset + "_#eta";

    ROOT::RDF::TH1DModel model_1D_eta(name_eta.c_str(), title_eta.c_str(), cfg.eta_plot.nbins, cfg.eta_plot.axis_min, 
    cfg.eta_plot.axis_max);

    std::string title_mll = ";" + cfg.mll_plot.title_axis + ";Efficiency;";
    std::string name_mll = cfg.general.dataset + "_m_{#mu+#mu-}";
    ROOT::RDF::TH1DModel model_1D_mll(name_mll.c_str(), title_mll.c_str(), cfg.mll_plot.nbins, cfg.mll_plot.axis_min, 
    cfg.mll_plot.axis_max);

    std::string title = ";" + cfg.eta_plot.title_axis + ";" + cfg.pt_plot.title_axis + ";";
    std::string name_plot = cfg.general.dataset + "_#eta VS p_{T}";
    ROOT::RDF::TH2DModel model_2D(name_plot.c_str(), title.c_str(), cfg.eta_plot.nbins, cfg.eta_plot.axis_min, cfg.eta_plot.axis_max, 
    cfg.pt_plot.nbins, cfg.pt_plot.axis_min, cfg.pt_plot.axis_max);
/*
    ROOT::RDF::TH2DModel model_2D_TP_Pt("h2_model_tp_pt", "; p_{T} Probe [GeV]; p_{T} Tag [GeV];", cfg.pt_plot.nbins, cfg.pt_plot.axis_min, cfg.pt_plot.axis_max,
    cfg.pt_plot.nbins, cfg.pt_plot.axis_min, cfg.pt_plot.axis_max);

    ROOT::RDF::TH2DModel model_2D_TP_Eta("h2_model_tp_eta", "; #eta Probe; #eta Tag;", cfg.eta_plot.nbins, cfg.eta_plot.axis_min, cfg.eta_plot.axis_max,
    cfg.eta_plot.nbins, cfg.eta_plot.axis_min, cfg.eta_plot.axis_max);

    ROOT::RDF::TH2DModel model_2D_RM_Pt("h2_model1", "; p_{T} gen [GeV]; p_{T} rec [GeV];", cfg.pt_plot.nbins, cfg.pt_plot.axis_min, cfg.pt_plot.axis_max,
    cfg.pt_plot.nbins, cfg.pt_plot.axis_min, cfg.pt_plot.axis_max);
    
    ROOT::RDF::TH2DModel model_2D_RM_Eta("h2_model2", "; #eta gen; #eta rec;", cfg.eta_plot.nbins, cfg.eta_plot.axis_min, cfg.eta_plot.axis_max,
    cfg.eta_plot.nbins, cfg.eta_plot.axis_min, cfg.eta_plot.axis_max);
*/
    std::vector<std::string> columns;

    if ((cfg.general.operation_mode.find("Selection") != std::string::npos) && (cfg.general.analysis_mode == "TagAndProbe")) {

        // ----------
        // Histograms
        // ----------

        auto h1_probe_pt = node.Histo1D(model_1D_pt, cfg.general.dataset + "_Probe_Pt");
        auto h1_probe_eta  = node.Histo1D(model_1D_eta, cfg.general.dataset + "_Probe_Eta");

        auto h1_mll = node.Histo1D(model_1D_mll, cfg.general.dataset + "_Mll");    

        auto h2_eta_pt = node.Histo2D(model_2D, cfg.general.dataset + "_Probe_Eta", cfg.general.dataset + "_Probe_Pt");

        // --------
        // Pipeline
        // --------

        AddToPipeline("Probe_Pt", h1_probe_pt);
        AddToPipeline("Probe_Eta", h1_probe_eta);
        AddToPipeline("InvMass", h1_mll);

        AddToPipeline("Eta vs Pt", h2_eta_pt);

        // -------------------
        // Saving Column names
        // -------------------

        std::vector<std::string> names = {cfg.general.dataset + "_Probe_Pt", cfg.general.dataset + "_Probe_Eta", 
            cfg.general.dataset + "_Mll", cfg.general.dataset + "_Mask_Pass"};

        columns.insert(columns.end(), names.begin(), names.end());

    } else if ((cfg.general.operation_mode.find("Selection") != std::string::npos) && 
    (cfg.general.analysis_mode == "ResponseMatrix")) {
        
        // ----------
        // Histograms
        // ----------
        
        auto h1_pt_gen = node.Histo1D(model_1D_pt, "Gen_Pt");
        auto h1_pt_rec = node.Histo1D(model_1D_pt, "Rec_Pt");

        auto h1_eta_gen = node.Histo1D(model_1D_eta, "Gen_Eta");
        auto h1_eta_rec = node.Histo1D(model_1D_eta, "Rec_Eta");

        //auto h2_pt_gen_rec = node.Histo2D(model_2D_RM_Pt, "Gen_Pt", "Rec_Pt");
        //auto h2_eta_gen_rec = node.Histo2D(model_2D_RM_Eta, "Gen_Eta", "Rec_Eta");

        // --------
        // Pipeline
        // --------

        AddToPipeline("P_{t} Generated distribution", h1_pt_gen);
        AddToPipeline("P_{t} Reconstructed distribution", h1_pt_rec);

        AddToPipeline("#eta Generated distribution", h1_eta_gen);
        AddToPipeline("#eta Reconstructed distribution", h1_eta_rec);

        //AddToPipeline("Response Matrix P_{t}", h2_pt_gen_rec);
        //AddToPipeline("Response Matrix #eta", h2_eta_gen_rec);
        
        std::vector<std::string> names = {"Gen_Pt", "Rec_Pt", "Gen_Eta", "Rec_Eta"};
        columns.insert(columns.end(), names.begin(), names.end());
    
    } else if ((cfg.general.operation_mode.find("Analysis") != std::string::npos) && (cfg.general.analysis_mode == "TagAndProbe")) {
        
        // Efficiency MonteCarlo
        std::vector<float> pt_bins = cfg.analysis.pt_bins;
        std::vector<float> eta_bins = cfg.analysis.eta_bins;

        ROOT::RDF::TH1DModel model_Eff_pt("h_Eff_pt", "Pt Efficiency; p_{T} [GeV]; Efficiency", pt_bins.size() - 1, pt_bins.data());
        ROOT::RDF::TH1DModel model_Eff_eta("h_Eff_eta", "Eta Efficiency; #eta; Efficiency", eta_bins.size() - 1, eta_bins.data());

        ROOT::RDF::TH2DModel model_2D("h2_Eff", "Efficiency map; #eta; p_{T} [GeV]", eta_bins.size() - 1, eta_bins.data(), 
        pt_bins.size() - 1, pt_bins.data());

        ROOT::RDF::RNode node_eff = node
            .Define(cfg.general.dataset + "_Probe_Pt_Pass", cfg.general.dataset + "_Probe_Pt[" + cfg.general.dataset + "_Mask_Pass]")
            .Define(cfg.general.dataset + "_Probe_Eta_Pass", cfg.general.dataset + "_Probe_Eta[" + cfg.general.dataset + "_Mask_Pass]");

        auto h1_probe_pt_pass = node_eff.Histo1D(model_Eff_pt, cfg.general.dataset + "_Probe_Pt_Pass");
        auto h1_eff_probe_pt_all = node_eff.Histo1D(model_Eff_pt, cfg.general.dataset + "_Probe_Pt");

        auto h1_probe_eta_pass = node_eff.Histo1D(model_Eff_eta, cfg.general.dataset + "_Probe_Eta_Pass");
        auto h1_eff_probe_eta_all = node_eff.Histo1D(model_Eff_eta, cfg.general.dataset + "_Probe_Eta");

        auto h2_probe_eta_pt_pass = node_eff.Histo2D(model_2D, cfg.general.dataset + "_Probe_Eta_Pass", cfg.general.dataset + "_Probe_Pt_Pass");
        auto h2_probe_eta_pt_all = node_eff.Histo2D(model_2D, cfg.general.dataset + "_Probe_Eta", cfg.general.dataset + "_Probe_Pt");

        AddToPipeline("Efficiency pt", h1_probe_pt_pass, h1_eff_probe_pt_all);
        AddToPipeline("Efficiency eta", h1_probe_eta_pass, h1_eff_probe_eta_all);
        AddToPipeline("Efficiency map", h2_probe_eta_pt_pass, h2_probe_eta_pt_all);
    }

    if (save_sel_data) {

        ROOT::RDF::RSnapshotOptions snapshot_opts;
        snapshot_opts.fMode = "UPDATE";
        snapshot_opts.fLazy = true;
        snapshot_opts.fOverwriteIfExists = true;

        std::string snaphot_name = cfg.general.dataset + "_" + cfg.general.analysis_mode + "_Tree";
        
        auto snapshot = node.Snapshot(snaphot_name, o_file_data, columns, snapshot_opts);
        
        std::cout << "Saving data selected from " << cfg.general.analysis_mode << " in file " << o_file_data << std::endl;
    
        snapshot_vec.push_back(snapshot);
    }
}


