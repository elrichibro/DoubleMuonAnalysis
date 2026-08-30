#include "Manager.h"

// #################
// ObjectTH1 Methods
// #################

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
    th1->Draw("AP");
}

// #################
// ObjectTH2 Methods
// #################

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

// #####################
// OutputManager Methods
// #####################

void OutputManager::AddToPipeline(const std::string& name, ROOT::RDF::RResultPtr<TH1D> hist) {
    pipeline.push_back(std::make_unique<ObjectTH1>(name, hist));
}

void OutputManager::AddToPipeline(const std::string& name, ROOT::RDF::RResultPtr<TH2D> hist) {
    pipeline.push_back(std::make_unique<ObjectTH2>(name, hist));
}

void OutputManager::Run() {
    std::unique_ptr<TFile> file = nullptr;
    
    if (save) {
        file = std::make_unique<TFile>(output_file.c_str(), "RECREATE");
    }

    for (auto& it : pipeline) {
        it->Process(); 

        // Saving
        if (save) {
            it->Write(*file); 
        }

        // Visualization
        if (visualize) {
            std::string c_name = "c_" + it->GetName();
            TCanvas* vis_canvas = new TCanvas(c_name.c_str(), it->GetName().c_str(), 800, 600);
            
            it->Draw(*vis_canvas);
            vis_canvas->Update();
        }
    }

    if (file && file->IsOpen()) {
        file->Close();
    }
}

void OutputManager::BookAnalysis(ROOT::RDF::RNode node, const config_struct& cfg, const std::string& mode) {
    
    std::string title_pt = ";" + cfg.pt_plot.title_axis + ";Efficiency;";

    ROOT::RDF::TH1DModel model_1D_pt("pt", title_pt.c_str(), cfg.pt_plot.nbins, cfg.pt_plot.axis_min, 
    cfg.pt_plot.axis_max);
    
    std::string title_eta = ";" + cfg.eta_plot.title_axis + ";Efficiency;";

    ROOT::RDF::TH1DModel model_1D_eta("eta", title_eta.c_str(), cfg.eta_plot.nbins, cfg.eta_plot.axis_min, 
    cfg.eta_plot.axis_max);

    std::string title = ";" + cfg.eta_plot.title_axis + ";" + cfg.pt_plot.title_axis + ";";

    ROOT::RDF::TH2DModel model_2D("h2", title.c_str(), cfg.eta_plot.nbins, cfg.eta_plot.axis_min, cfg.eta_plot.axis_max, 
    cfg.pt_plot.nbins, cfg.pt_plot.axis_min, cfg.pt_plot.axis_max);
    
    if (mode == "TagAndProbe") {
        auto h_pt_num = node.Histo1D(model_1D_pt, "Probe_Pt_Pass");
        auto h_pt_den  = node.Histo1D(model_1D_pt, "Probe_Pt_All");

        auto h_eta_num = node.Histo1D(model_1D_eta, "Probe_Eta_Pass");
        auto h_eta_den  = node.Histo1D(model_1D_eta, "Probe_Eta_All");

        auto h2_den = node.Histo2D(model_2D, "Probe_Eta_All", "Probe_Pt_All");
        auto h2_num = node.Histo2D(model_2D, "Probe_Eta_Pass", "Probe_Pt_Pass");

        //AddToPipeline("Pt_Pass", h_pt_num);
        //AddToPipeline("Pt_Total", h_pt_den);

        //AddToPipeline("Eta_Pass", h_eta_num);
        //AddToPipeline("Eta_Total", h_eta_den);
        
        AddToPipeline("Efficiency pt", h_pt_num, h_pt_den);
        AddToPipeline("Efficiency eta", h_eta_num, h_eta_den);
        AddToPipeline("Efficiency map", h2_num, h2_den);
    
    } else if (mode == "ResponseMatrix") {

    }
}


