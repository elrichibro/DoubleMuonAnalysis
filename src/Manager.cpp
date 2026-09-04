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
    
    std::unique_ptr<TFile> file_plots = nullptr;
    if (save_sel_plots) {
        file_plots = std::make_unique<TFile>(o_file_plots.c_str(), "RECREATE");
        std::cout << "Save file booked." << std::endl;
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

    for (auto& snap : snapshot_vec) {
        snap.GetValue(); 
        
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> DeltaT = end_time - start_time;
        
        std::cout << "Snapshot time: " << DeltaT.count() << std::endl;
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
    std::string title_eta = ";" + cfg.eta_plot.title_axis + ";Efficiency;";
    std::string title_mll = ";" + cfg.mll_plot.title_axis + ";Efficiency;";
    std::string title = ";" + cfg.eta_plot.title_axis + ";" + cfg.pt_plot.title_axis + ";";

    ROOT::RDF::TH1DModel model_1D_pt("p_{T} stats", title_pt.c_str(), cfg.pt_plot.nbins, cfg.pt_plot.axis_min, 
    cfg.pt_plot.axis_max);

    ROOT::RDF::TH1DModel model_1D_eta("#eta stats", title_eta.c_str(), cfg.eta_plot.nbins, cfg.eta_plot.axis_min, 
    cfg.eta_plot.axis_max);

    ROOT::RDF::TH1DModel model_1D_mll("m_{#mu+#mu-} stats", title_mll.c_str(), cfg.mll_plot.nbins, cfg.mll_plot.axis_min, 
    cfg.mll_plot.axis_max);

    ROOT::RDF::TH2DModel model_2D("Stats", title.c_str(), cfg.eta_plot.nbins, cfg.eta_plot.axis_min, cfg.eta_plot.axis_max, 
    cfg.pt_plot.nbins, cfg.pt_plot.axis_min, cfg.pt_plot.axis_max);

    ROOT::RDF::TH2DModel model_2D_TP_Pt("h2_model_tp_pt", "; p_{T} Probe [GeV]; p_{T} Tag [GeV];", cfg.pt_plot.nbins, cfg.pt_plot.axis_min, cfg.pt_plot.axis_max,
    cfg.pt_plot.nbins, cfg.pt_plot.axis_min, cfg.pt_plot.axis_max);

    ROOT::RDF::TH2DModel model_2D_TP_Eta("h2_model_tp_eta", "; #eta Probe; #eta Tag;", cfg.eta_plot.nbins, cfg.eta_plot.axis_min, cfg.eta_plot.axis_max,
    cfg.eta_plot.nbins, cfg.eta_plot.axis_min, cfg.eta_plot.axis_max);

    ROOT::RDF::TH2DModel model_2D_RM_Pt("h2_model1", "; p_{T} gen [GeV]; p_{T} rec [GeV];", cfg.pt_plot.nbins, cfg.pt_plot.axis_min, cfg.pt_plot.axis_max,
    cfg.pt_plot.nbins, cfg.pt_plot.axis_min, cfg.pt_plot.axis_max);
    
    ROOT::RDF::TH2DModel model_2D_RM_Eta("h2_model2", "; #eta gen; #eta rec;", cfg.eta_plot.nbins, cfg.eta_plot.axis_min, cfg.eta_plot.axis_max,
    cfg.eta_plot.nbins, cfg.eta_plot.axis_min, cfg.eta_plot.axis_max);

    std::vector<std::string> columns;

    if ((cfg.general.operation_mode.find("Selection") != std::string::npos) && (cfg.general.analysis_mode == "TagAndProbe")) {

        // For TEfficiency -> just visualization -> NOT analysis
        ROOT::RDF::RNode node_eff = node
            .Define("Probe_Pt_All", [](const ROOT::RVec<float>& pass, const ROOT::RVec<float>& fail) {
                return ROOT::VecOps::Concatenate(pass, fail);
            }, {"Probe_Pt_Pass", "Probe_Pt_Fail"})
            .Define("Probe_Eta_All", [](const ROOT::RVec<float>& pass, const ROOT::RVec<float>& fail) {
                return ROOT::VecOps::Concatenate(pass, fail);
            }, {"Probe_Eta_Pass", "Probe_Eta_Fail"});

        // ----------
        // Histograms
        // ----------

        
        // Pt/Eta/Mll Pass/Fail
        /*
        auto h1_probe_pt_pass = node_eff.Histo1D(model_1D_pt, "Probe_Pt_Pass");
        auto h1_probe_pt_fail = node_eff.Histo1D(model_1D_pt, "Probe_Pt_Fail");

        auto h1_probe_eta_pass = node_eff.Histo1D(model_1D_eta, "Probe_Eta_Pass");
        auto h1_probe_eta_fail  = node_eff.Histo1D(model_1D_eta, "Probe_Eta_Fail");

        auto h1_mll_pass = node_eff.Histo1D(model_1D_mll, "Mll_Pass");
        auto h1_mll_fail  = node_eff.Histo1D(model_1D_mll, "Mll_Fail");
        
        // Tag
        auto h1_tag_pt_pass = node_eff.Histo1D(model_1D_pt, "Tag_Pt_Pass");
        auto h1_tag_eta_pass = node_eff.Histo1D(model_1D_eta, "Tag_Eta_Pass");

        // TEff
        auto h1_eff_probe_pt_all = node_eff.Histo1D(model_1D_pt, "Probe_Pt_All");
        auto h1_eff_probe_eta_all = node_eff.Histo1D(model_1D_eta, "Probe_Eta_All");

        auto h2_eta_pt_pass = node_eff.Histo2D(model_2D, "Probe_Eta_Pass", "Probe_Pt_Pass");
        auto h2_eta_pt_fail = node_eff.Histo2D(model_2D, "Probe_Eta_Fail", "Probe_Pt_Fail");
        auto h2_eff_eta_pt_all = node_eff.Histo2D(model_2D, "Probe_Eta_All", "Probe_Pt_All");        

        auto h2_probe_tag_pt_fail = node_eff.Histo2D(model_2D_TP_Pt, "Probe_Pt_Fail", "Tag_Pt_Fail");
        auto h2_probe_tag_pt_pass = node_eff.Histo2D(model_2D_TP_Pt, "Probe_Pt_Pass", "Tag_Pt_Pass");

        auto h2_probe_tag_eta_fail = node_eff.Histo2D(model_2D_TP_Eta, "Probe_Eta_Fail", "Tag_Eta_Fail");
        auto h2_probe_tag_eta_pass = node_eff.Histo2D(model_2D_TP_Eta, "Probe_Eta_Pass", "Tag_Eta_Pass");
        */

        // --------
        // Pipeline
        // --------

        /*
        // Probes
        AddToPipeline("Probe_Pt_Pass", h1_probe_pt_pass);
        AddToPipeline("Probe_Pt_Fail", h1_probe_pt_fail);

        AddToPipeline("Probe_Eta_Pass", h1_probe_eta_pass);
        AddToPipeline("Probe_Eta_Fail", h1_probe_eta_fail);

        AddToPipeline("InvMass_Pass", h1_mll_pass);
        AddToPipeline("InvMass_Fail", h1_mll_fail);

        // Tag
        AddToPipeline("Tag_Pt_Pass", h1_tag_pt_pass);
        AddToPipeline("Tag_Eta_Pass", h1_tag_eta_pass);

        AddToPipeline("Efficiency pt", h1_probe_pt_pass, h1_eff_probe_pt_all);
        AddToPipeline("Efficiency eta", h1_probe_eta_pass, h1_eff_probe_eta_all);

        AddToPipeline("Correlaton Tag/Probe Pt Pass", h2_probe_tag_pt_pass);
        AddToPipeline("Correlaton Tag/Probe Pt Fail", h2_probe_tag_pt_fail);
        AddToPipeline("Correlaton Tag/Probe Eta Pass", h2_probe_tag_eta_pass);
        AddToPipeline("Correlaton Tag/Probe Eta Fail", h2_probe_tag_eta_fail);

        AddToPipeline("Efficiency map", h2_eta_pt_pass, h2_eff_eta_pt_all);
        */

        std::vector<std::string> names = {"Probe_Pt_Pass", "Probe_Pt_Fail", "Probe_Eta_Pass", 
        "Probe_Eta_Fail", "Mll_Pass", "Mll_Fail"};

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

        auto h2_pt_gen_rec = node.Histo2D(model_2D_RM_Pt, "Gen_Pt", "Rec_Pt");
        auto h2_eta_gen_rec = node.Histo2D(model_2D_RM_Eta, "Gen_Eta", "Rec_Eta");

        // --------
        // Pipeline
        // --------

        AddToPipeline("P_{t} Generated distribution", h1_pt_gen);
        AddToPipeline("P_{t} Reconstructed distribution", h1_pt_rec);

        AddToPipeline("#eta Generated distribution", h1_eta_gen);
        AddToPipeline("#eta Reconstructed distribution", h1_eta_rec);

        AddToPipeline("Response Matrix P_{t}", h2_pt_gen_rec);
        AddToPipeline("Response Matrix #eta", h2_eta_gen_rec);
        
        std::vector<std::string> names = {"Gen_Pt", "Rec_Pt", "Gen_Eta", "Rec_Eta"};
        columns.insert(columns.end(), names.begin(), names.end());
    
    } else if ((cfg.general.operation_mode.find("Analysis") != std::string::npos) && 
    (cfg.general.analysis_mode == "TagAndProbe")) {
        
        // Efficiency MonteCarlo
        std::vector<float> pt_bins = cfg.analysis.pt_bins;
        std::vector<float> eta_bins = cfg.analysis.eta_bins;

        ROOT::RDF::TH1DModel model_Eff_pt("h_Eff_pt", "Pt Efficiency; p_{T} [GeV]; Efficiency", pt_bins.size() - 1, pt_bins.data());
        ROOT::RDF::TH1DModel model_Eff_eta("h_Eff_eta", "Eta Efficiency; #eta; Efficiency", eta_bins.size() - 1, eta_bins.data());

        ROOT::RDF::RNode node_eff = node
            .Define("Probe_Pt_All", [](const ROOT::RVec<float>& pass, const ROOT::RVec<float>& fail) {
                return ROOT::VecOps::Concatenate(pass, fail);
            }, {"Probe_Pt_Pass", "Probe_Pt_Fail"})
            .Define("Probe_Eta_All", [](const ROOT::RVec<float>& pass, const ROOT::RVec<float>& fail) {
                return ROOT::VecOps::Concatenate(pass, fail);
            }, {"Probe_Eta_Pass", "Probe_Eta_Fail"});

        auto h1_probe_pt_pass = node_eff.Histo1D(model_Eff_pt, "Probe_Pt_Pass");
        auto h1_eff_probe_pt_all = node_eff.Histo1D(model_Eff_pt, "Probe_Pt_All");

        auto h1_probe_eta_pass = node_eff.Histo1D(model_Eff_eta, "Probe_Eta_Pass");
        auto h1_eff_probe_eta_all = node_eff.Histo1D(model_Eff_eta, "Probe_Eta_All");

        AddToPipeline("Efficiency pt", h1_probe_pt_pass, h1_eff_probe_pt_all);
        AddToPipeline("Efficiency eta", h1_probe_eta_pass, h1_eff_probe_eta_all);
    }

    if (save_sel_data) {
        ROOT::RDF::RSnapshotOptions snapshot_opts;
        snapshot_opts.fMode = "UPDATE";
        snapshot_opts.fLazy = true;
        snapshot_opts.fOverwriteIfExists = true;

        std::string snaphot_name = cfg.general.dataset + cfg.general.analysis_mode + "_Tree";
        
        auto snapshot = node.Snapshot(snaphot_name, o_file_data, columns, snapshot_opts);
        
        std::cout << "Saving data selected from " << cfg.general.analysis_mode << " in file " << o_file_data << std::endl;
    
        snapshot_vec.push_back(snapshot);
    }
}


