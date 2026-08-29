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

void ObjectTEff::Process(TCanvas& canvas, TFile& file) {
    auto eff_pt = std::make_unique<TEfficiency>(*hist_num, *hist_den);
    eff_pt->Write();
    eff_pt->Draw();
    canvas.SaveAs((name_eff + ".root").c_str());
}

void OutputManager::AddToPipeline(const std::string& name, ROOT::RDF::RResultPtr<TH1D> hist) {
    pipeline.push_back(std::make_unique<ObjectTH1>(name, hist));
}

void OutputManager::AddToPipeline(const std::string& name, ROOT::RDF::RResultPtr<TH2D> hist) {
    pipeline.push_back(std::make_unique<ObjectTH2>(name, hist));
}

void OutputManager::AddToPipeline(const std::string& name, ROOT::RDF::RResultPtr<TH1D> hist1, ROOT::RDF::RResultPtr<TH1D> hist2) {
    pipeline.push_back(std::make_unique<ObjectTEff>(name, hist1, hist2));
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