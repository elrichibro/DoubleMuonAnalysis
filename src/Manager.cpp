#include "Manager.h"

#include <string>
#include <vector>
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
    for (auto& snap : snapshot_vec) {
        snap.GetValue(); 
    }
    
    
    std::unique_ptr<TFile> file_plots = nullptr;
    if (save_sel_plots) {
        file_plots = std::make_unique<TFile>(o_file_plots.c_str(), "RECREATE");
        std::cout << "Save file booked." << std::endl;
    }

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
    }

    if (file_plots && file_plots->IsOpen()) {
        file_plots->Close();
    }
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

    ROOT::RDF::TH2DModel model_2D_1("h2_model1", "; p_{T} gen [GeV]; p_{T} rec [GeV];", cfg.pt_plot.nbins, cfg.pt_plot.axis_min, cfg.pt_plot.axis_max,
    cfg.pt_plot.nbins, cfg.pt_plot.axis_min, cfg.pt_plot.axis_max);
    
    ROOT::RDF::TH2DModel model_2D_2("h2_model2", "; #eta gen; #eta rec;", cfg.eta_plot.nbins, cfg.eta_plot.axis_min, cfg.eta_plot.axis_max,
    cfg.eta_plot.nbins, cfg.eta_plot.axis_min, cfg.eta_plot.axis_max);

    std::vector<std::string> columns;

    if (cfg.general.analysis_mode == "TagAndProbe") {

        // ----------
        // Histograms
        // ----------

        auto h1_pt_num = node.Histo1D(model_1D_pt, "Probe_Pt_Pass");
        auto h1_pt_den  = node.Histo1D(model_1D_pt, "Probe_Pt_All");

        auto h1_eta_num = node.Histo1D(model_1D_eta, "Probe_Eta_Pass");
        auto h1_eta_den  = node.Histo1D(model_1D_eta, "Probe_Eta_All");

        auto h1_mll_num = node.Histo1D(model_1D_mll, "Probe_Mll_Pass");
        auto h1_mll_den  = node.Histo1D(model_1D_mll, "Probe_Mll_All");

        auto h2_eta_pt_den = node.Histo2D(model_2D, "Probe_Eta_All", "Probe_Pt_All");
        auto h2_eta_pt_num = node.Histo2D(model_2D, "Probe_Eta_Pass", "Probe_Pt_Pass");

        // --------
        // Pipeline
        // --------

        AddToPipeline("Pt_Pass", h1_pt_num);
        AddToPipeline("Pt_Total", h1_pt_den);

        //AddToPipeline("Eta_Pass", h1_eta_num);
        //AddToPipeline("Eta_Total", h1_eta_den);

        AddToPipeline("InvMass_Total", h1_mll_den);
        AddToPipeline("InvMass_Pass", h1_mll_num);

        AddToPipeline("Efficiency pt", h1_pt_num, h1_pt_den);
        AddToPipeline("Efficiency eta", h1_eta_num, h1_eta_den);
        
        AddToPipeline("Efficiency map", h2_eta_pt_num, h2_eta_pt_den);
    
        std::vector<std::string> names = {"Probe_Pt_Pass", "Probe_Pt_All", "Probe_Eta_Pass", 
        "Probe_Eta_All", "Probe_Mll_Pass", "Probe_Mll_All"};
        
        columns.insert(columns.end(), names.begin(), names.end());

    } else if (cfg.general.analysis_mode == "ResponseMatrix") {
        
        // ----------
        // Histograms
        // ----------
        
        auto h1_pt_gen = node.Histo1D(model_1D_pt, "Gen_Pt");
        auto h1_pt_rec = node.Histo1D(model_1D_pt, "Rec_Pt");

        auto h1_eta_gen = node.Histo1D(model_1D_eta, "Gen_Eta");
        auto h1_eta_rec = node.Histo1D(model_1D_eta, "Rec_Eta");

        auto h2_pt_gen_rec = node.Histo2D(model_2D_1, "Gen_Pt", "Rec_Pt");
        auto h2_eta_gen_rec = node.Histo2D(model_2D_2, "Gen_Eta", "Rec_Eta");

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
    }

    if (save_sel_data) {
        ROOT::RDF::RSnapshotOptions snapshot_opts;
        snapshot_opts.fMode = "UPDATE";
        snapshot_opts.fLazy = true;
        snapshot_opts.fOverwriteIfExists = true;

        auto snapshot = node.Snapshot(cfg.general.analysis_mode + "_Tree", o_file_data, columns, snapshot_opts);
        std::cout << "Saving data selected from " << cfg.general.analysis_mode << " in file " << o_file_data << std::endl;
    
        snapshot_vec.push_back(snapshot);
    }
}


