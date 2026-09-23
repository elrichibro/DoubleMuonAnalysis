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

void OutputSelManager::BookAction(ROOT::RDF::RNode node, const config_struct& cfg) {
    
    if ((cfg.general.operation_mode.find("Selection") != std::string::npos) && (cfg.selection.selection_mode == "RespMatrix")) {
        // -------------------
        // Saving Column names
        // -------------------

        std::vector<std::string> names = {"Matched", "Missed", "Faked", "Rec_InvMass", "Gen_InvMass", "Rec_Pt", "Gen_Pt", "Rec_Y", "Gen_Y", "Rec_Phis", "Gen_Phis"};
        column_names.insert(column_names.end(), names.begin(), names.end());
    
    } else if ((cfg.general.operation_mode.find("Selection") != std::string::npos) && (cfg.selection.selection_mode == "Event")) {
        // -------------------
        // Saving Column names
        // -------------------

        std::vector<std::string> names = {"InvariantMass", "Pt_Z", "Y_Z", "Phis_Z"};
        column_names.insert(column_names.end(), names.begin(), names.end());
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