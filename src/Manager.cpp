#include "Manager.h"

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
    th2->Draw("AP");
}


void OutputManager::AddToPipeline(const std::string& name, ROOT::RDF::RResultPtr<TH1D> hist) {
    pipeline.push_back(std::make_unique<ObjectTH1>(name, hist));
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

