#include "Manager.h"

void ObjectTH1::Process(TCanvas& canvas, TFile& file) {
    th1->Write();
    th1->Draw();
    canvas.SaveAs((name_th1 + ".root").c_str());
}

void ObjectTH2::Process(TCanvas& canvas, TFile& file) {
    th2->Write();
    th2->Draw();
    canvas.SaveAs((name_th2 + ".root").c_str());
}

void OutputManager::AddToPipeline(const std::string& name, ROOT::RDF::RResultPtr<TH1D> hist) {
    pipeline.push_back(std::make_unique<ObjectTH1>(name, hist));
}

void OutputManager::Run(const std::string& filename) {
    TFile file(filename.c_str(), "RECREATE");
    TCanvas canvas("c", "Canvas", 800, 600);
    
    for (auto& it : pipeline) {
        canvas.Clear();
        it->Process(canvas, file);
    }
    
    file.Close();
}