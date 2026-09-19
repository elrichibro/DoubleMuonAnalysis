#include "Manager.h"

#include <string>
#include <vector>

#include <chrono>

#include "Unfold.h"
/*
Classes:
    - ObjectTH1
    - ObjectTH2
    - ObjectTEff

    - OutputManager
*/

// ------------------------------------------------------------------------------------------------------------------------------------

void ObjectTH1::Write(TFile& file) {
    if (th1) {
        file.cd();
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

void ObjectTH2::Write(TFile& file) {
    if (th2) {
        file.cd();
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

void OutputSelManager::AddToPipeline(const std::string& name, ROOT::RDF::RResultPtr<TH1D> hist) {
    pipeline.push_back(std::make_unique<ObjectTH1>(name, hist));
}

void OutputSelManager::AddToPipeline(const std::string& name, ROOT::RDF::RResultPtr<TH2D> hist) {
    pipeline.push_back(std::make_unique<ObjectTH2>(name, hist));
}

void OutputSelManager::Run() {
    // Check execution time - START
    auto start_time = std::chrono::high_resolution_clock::now();

    for (auto& snap : snapshot_vec) {
        // if snapshot is enabled starts the Event Loop
        snap.GetValue(); 
        
        // Check execution time - STOP
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> DeltaT = end_time - start_time;
        
        std::cout << "Snapshot time: " << DeltaT.count() << std::endl;
    }

    std::unique_ptr<TFile> file_plots = nullptr;
    
    if ((config.selection.save_sel_plots) && (config.general.operation_mode == "Selection")) {
        file_plots = std::make_unique<TFile>(config.selection.o_sel_file_plots.c_str(), "UPDATE");
        std::cout << "Saving plots in file: " <<  config.selection.o_sel_file_plots << std::endl;
    
    } else if ((config.selection.save_sel_plots) && (config.general.operation_mode == "Analysis")) {
        file_plots = std::make_unique<TFile>(config.analysis.o_fit_file.c_str(), "UPDATE");
        std::cout << "Saving plots in file: " <<  config.analysis.o_fit_file << std::endl;
    
    }

    int i = 0;

    for (auto& it : pipeline) {
        // if the First PipelineObj is a TEfficiency obj -> starts the Event Loop
        it->Process();

        // Saving plots
        if (config.selection.save_sel_plots) {
            // if save plots is enabled starts the Event Loop
            it->Write(*file_plots);
        }

        // Visualization option
        if (config.selection.visual_sel) {

            std::string c_name = "c_" + it->GetName();
            TCanvas* vis_canvas = new TCanvas(c_name.c_str(), it->GetName().c_str(), config.canvas.width, config.canvas.height);
            
            it->Draw(*vis_canvas);
            vis_canvas->Update();
        }

        // Check execution time - STOP
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> DeltaT = end_time - start_time;

        std::cout << "PipelineObj number: " << i << ", execution time: " << DeltaT.count() << std::endl;
        i++;
    }

    if (file_plots && file_plots->IsOpen()) {
        file_plots->Close();
        
        // Check execution time - STOP
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> DeltaT = end_time - start_time;
        
        std::cout << "Output plots file closed at: " << DeltaT.count() << std::endl;
    }

    std::cout << "Ending of OutputManager::Run()" << std::endl;
}

// ------------------------------------------------------------------------------------------------------------------------------------

void OutputSelManager::BookAnalysis(ROOT::RDF::RNode node, const config_struct& cfg) {
    
    // ---------------------------------
    // Defining general histogram models
    // ---------------------------------

    ROOT::RDF::TH2DModel model_2D_TP_Pt("h2_model_tp_pt", "; p_{T} Probe [GeV]; p_{T} Tag [GeV];", cfg.pt_plot.nbins, cfg.pt_plot.axis_min, cfg.pt_plot.axis_max,
    cfg.pt_plot.nbins, cfg.pt_plot.axis_min, cfg.pt_plot.axis_max);

    ROOT::RDF::TH2DModel model_2D_TP_Eta("h2_model_tp_eta", "; #eta Probe; #eta Tag;", cfg.eta_plot.nbins, cfg.eta_plot.axis_min, cfg.eta_plot.axis_max,
    cfg.eta_plot.nbins, cfg.eta_plot.axis_min, cfg.eta_plot.axis_max);

    ROOT::RDF::TH2DModel model_2D_RM_Pt("h2_model1", "; p_{T} gen [GeV]; p_{T} rec [GeV];", cfg.pt_plot.nbins, cfg.pt_plot.axis_min, cfg.pt_plot.axis_max,
    cfg.pt_plot.nbins, cfg.pt_plot.axis_min, cfg.pt_plot.axis_max);
    
    ROOT::RDF::TH2DModel model_2D_RM_Eta("h2_model2", "; #eta gen; #eta rec;", cfg.eta_plot.nbins, cfg.eta_plot.axis_min, cfg.eta_plot.axis_max,
    cfg.eta_plot.nbins, cfg.eta_plot.axis_min, cfg.eta_plot.axis_max);

    // vector needed for Snapshot operation

    if ((cfg.general.operation_mode.find("Selection") != std::string::npos) && (cfg.selection.selection_mode == "TagAndProbe")) {
        
        // ------------
        // Histo Models
        // ------------

        std::string title_pt = ";" + cfg.pt_plot.title_axis + ";Efficiency;";
        std::string name_pt = cfg.selection.dataset + "_p_{T}";
        ROOT::RDF::TH1DModel model_1D_pt(name_pt.c_str(), title_pt.c_str(), cfg.pt_plot.nbins, cfg.pt_plot.axis_min, 
        cfg.pt_plot.axis_max);

        std::string title_eta = ";" + cfg.eta_plot.title_axis + ";Efficiency;";
        std::string name_eta = cfg.selection.dataset + "_#eta";
        ROOT::RDF::TH1DModel model_1D_eta(name_eta.c_str(), title_eta.c_str(), cfg.eta_plot.nbins, cfg.eta_plot.axis_min, 
        cfg.eta_plot.axis_max);

        std::string title_mll = ";" + cfg.mll_plot.title_axis + ";Efficiency;";
        std::string name_mll = cfg.selection.dataset + "_m_{#mu+#mu-}";
        ROOT::RDF::TH1DModel model_1D_mll(name_mll.c_str(), title_mll.c_str(), cfg.mll_plot.nbins, cfg.mll_plot.axis_min, 
        cfg.mll_plot.axis_max);

        std::string title = ";" + cfg.eta_plot.title_axis + ";" + cfg.pt_plot.title_axis + ";";
        std::string name_plot = cfg.selection.dataset + "_#eta VS p_{T}";
        ROOT::RDF::TH2DModel model_2D(name_plot.c_str(), title.c_str(), cfg.eta_plot.nbins, cfg.eta_plot.axis_min, cfg.eta_plot.axis_max, 
        cfg.pt_plot.nbins, cfg.pt_plot.axis_min, cfg.pt_plot.axis_max);

        // ----------
        // Histograms
        // ----------

        auto h1_probe_pt = node.Histo1D(model_1D_pt, cfg.selection.dataset + "_Probe_Pt");
        auto h1_probe_eta  = node.Histo1D(model_1D_eta, cfg.selection.dataset + "_Probe_Eta");

        auto h1_mll = node.Histo1D(model_1D_mll, cfg.selection.dataset + "_Mll");    

        auto h2_eta_pt = node.Histo2D(model_2D, cfg.selection.dataset + "_Probe_Eta", cfg.selection.dataset + "_Probe_Pt");

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

        std::vector<std::string> names = {cfg.selection.dataset + "_Probe_Pt", cfg.selection.dataset + "_Probe_Eta", 
            cfg.selection.dataset + "_Mll", cfg.selection.dataset + "_Mask_Pass"};

        column_names.insert(column_names.end(), names.begin(), names.end());

    } else if ((cfg.general.operation_mode.find("Selection") != std::string::npos) && (cfg.selection.selection_mode == "RespMatrix")) {
        
        std::string title_pt = ";" + cfg.pt_plot.title_axis + ";Efficiency;";
        std::string name_pt = cfg.selection.dataset + "_p_{T}";
        ROOT::RDF::TH1DModel model_1D_pt(name_pt.c_str(), title_pt.c_str(), cfg.pt_plot.nbins, cfg.pt_plot.axis_min, 
        cfg.pt_plot.axis_max);

        RespMatrixHisto histo = BuildRespMatrixHisto(node, cfg);

        // --------
        // Pipeline
        // --------

        AddToPipeline("P_{t, Z0} Response Matrix", histo.histo_pt);
        AddToPipeline("Y_{Z0}", histo.histo_y);
        AddToPipeline("#Phi_{Z0}^{*}", histo.histo_phis);
 
        // -------------------
        // Saving Column names
        // -------------------

        std::vector<std::string> names = {"Rec_InvMass", "Gen_InvMass", "Rec_Pt", "Gen_Pt", "Rec_Y", "Gen_Y", "Rec_Phis", "Gen_Phis"};
        column_names.insert(column_names.end(), names.begin(), names.end());
    
    } else if ((cfg.general.operation_mode.find("Analysis") != std::string::npos) && (cfg.analysis.analysis_mode == "TagAndProbe_MC")) {
        
        // Efficiency MonteCarlo
        std::vector<float> pt_bins = cfg.templ.pt_bins;
        std::vector<float> eta_bins = cfg.templ.eta_bins;

        ROOT::RDF::TH1DModel model_Eff_pt("h_Eff_pt", "Pt Efficiency; p_{T} [GeV]; Efficiency", pt_bins.size() - 1, pt_bins.data());
        ROOT::RDF::TH1DModel model_Eff_eta("h_Eff_eta", "Eta Efficiency; #eta; Efficiency", eta_bins.size() - 1, eta_bins.data());

        ROOT::RDF::TH2DModel model_2D("h2_Eff", "Efficiency map; #eta; p_{T} [GeV]", eta_bins.size() - 1, eta_bins.data(), 
        pt_bins.size() - 1, pt_bins.data());

        ROOT::RDF::RNode node_eff = node
            .Define("MC_Probe_Pt_Pass", "MC_Probe_Pt[MC_Mask_Pass]")
            .Define("MC_Probe_Eta_Pass", "MC_Probe_Eta[MC_Mask_Pass]");

        auto h1_probe_pt_pass = node_eff.Histo1D(model_Eff_pt, "MC_Probe_Pt_Pass");
        auto h1_eff_probe_pt_all = node_eff.Histo1D(model_Eff_pt, "MC_Probe_Pt");

        auto h1_probe_eta_pass = node_eff.Histo1D(model_Eff_eta, "MC_Probe_Eta_Pass");
        auto h1_eff_probe_eta_all = node_eff.Histo1D(model_Eff_eta, "MC_Probe_Eta");

        auto h2_probe_eta_pt_pass = node_eff.Histo2D(model_2D, "MC_Probe_Eta_Pass", "MC_Probe_Pt_Pass");
        auto h2_probe_eta_pt_all = node_eff.Histo2D(model_2D, "MC_Probe_Eta", "MC_Probe_Pt");

        AddToPipeline("Efficiency pt", h1_probe_pt_pass, h1_eff_probe_pt_all);
        AddToPipeline("Efficiency eta", h1_probe_eta_pass, h1_eff_probe_eta_all);
        AddToPipeline("Efficiency map", h2_probe_eta_pt_pass, h2_probe_eta_pt_all);
    }

    // ---------------
    // Snapshot option
    // ---------------

    if (config.selection.save_sel_data) {

        ROOT::RDF::RSnapshotOptions snapshot_opts;
        snapshot_opts.fMode = "UPDATE";
        snapshot_opts.fLazy = true;
        snapshot_opts.fOverwriteIfExists = true;

        std::string snaphot_name = cfg.selection.dataset + "_" + cfg.selection.selection_mode + "_Tree";
        
        auto snapshot = node.Snapshot(snaphot_name, config.selection.o_sel_file_data, column_names, snapshot_opts);
        
        std::cout << "Saving data selected from " << cfg.selection.selection_mode << " in file " << config.selection.o_sel_file_data << std::endl;
    
        snapshot_vec.push_back(snapshot);// Needed for scope visibility -> Smart pointer for Event Loop action
    }
}