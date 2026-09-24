#include <ROOT/RDataFrame.hxx>

#include <iostream>
#include <string>

#include <TApplication.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMarker.h>
#include "TEfficiency.h"

#include "Config.h"
#include "Utils.h"
#include "Manager.h"
#include "Unfold.h"
#include "Selection.h"
#include "Event.h"
#include "ControlCheck.h"

int main(int argc, char* argv[]) {

    config_struct cfg;
    std::string json_path = "";
    
    int verbose = 0;
    bool control = false;
    bool visualize = false;
    
    bool save_plots = false;
    bool save_data = false;
    
    // ---------
    // INTERFACE
    // ---------
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg.rfind(".json") != std::string::npos) {
            json_path = arg;
        } else if (arg == "--verbose" || arg == "-v") {
            verbose = 1;
        } else if ((arg == "--visualize") || (arg == "-vis")) {
            visualize = true;
        } else if ((arg == "--control") || (arg == "-c")) {
            control = true;
        } else {
            std::cout << "ERROR: invalid input command, please try again, exinting.\n" << std::endl;
            return 1;
        }
    }

    // ------------------------------------------------------------------------------------------------------------------------------------
    // JSON CONFIGURATION
    // ------------------------------------------------------------------------------------------------------------------------------------

    if (Configure(cfg, json_path) != 0) {
        std::cout << "Configurations fails, check needed, exiting." << std::endl;
        return 1;
    }

    verbose = cfg.general.verbose;
    visualize = cfg.general.visualize;

    // Verbose JSON configuration
    if (verbose) {
        Verbose_config(cfg);
    }

    if (control) {
        std::cout << "Control config settup done, exiting." << std::endl;
        return 0;
    }

    const flags_config flags_RM = cfg.flag_RM;
    const cuts_config cuts_RM = cfg.cut_RM;

    validation_type validation_map = Validation_load(cfg.io.val_file);
    std::cout << "Validation Map created." << std::endl;

    std::string dataset_tree = "";
    std::string dataset_file = "";

    ROOT::EnableImplicitMT();// MultiThread option: ON

    // ------------------------------------------------------------------------------------------------------------------------------------
    // VISUALIZATION OPTION
    // ------------------------------------------------------------------------------------------------------------------------------------
    
    TApplication* app = nullptr;
    if (visualize) {
        std::cout << "TApplication inializating ..." << std::endl;
        app = new TApplication("app", &argc, argv);
    }

    // ----------
    // Acceptance
    // ----------

    if (cfg.general.operation_mode.find("Acceptance") != std::string::npos) {

        if (cfg.acceptance.dataset == "DATA") {
            std::cout << "ERROR: Invalid dataset with Acceptance Operation Mode, only MC is permitted, exiting..." << std::endl;
            return 1;
        } else if (cfg.acceptance.dataset == "MC") {
            dataset_tree = cfg.io.tree_mc_name;
            dataset_file = cfg.io.in_mc_file;
        }

        ROOT::EnableImplicitMT();// MultiThread option: ON

        ROOT::RDataFrame data_frame(dataset_tree, dataset_file);

        if (verbose){ 
            std::cout << "RDataFrame object created, unpacking tree: " << dataset_tree 
            << ", from file: " << dataset_file << ", starting selection ..." << std::endl;
        }
            
        ROOT::RDF::RNode node_ACC = data_frame;

        std::vector<float> results = CalculateAcceptance(node_ACC, "aFSR", 2);
        
        if (results.size() != 2) {
            std::cout << "ERROR: invalid results size: " << results.size() << ", exiting.." << std::endl;
            return 1;
        }
            
        std::cout << "Geometrical acceptance: " << results.at(0) << "+-" << results.at(1) << std::endl;
    }

    // ----------
    // Resolution
    // ----------

    if (cfg.general.operation_mode.find("Resolution") != std::string::npos) {

        ROOT::RDataFrame data_frame("MC_RespMatrix_Tree", cfg.selection.o_sel_file_data);
        ROOT::RDF::RNode node_resolution = data_frame;

        //ROOT::EnableImplicitMT();// MultiThread option: ON

        ResolutionResults resolution =  CalculateResolution(node_resolution, cfg);
        
        if (app != nullptr) {                
            std::cout << "Starting visualization..." << std::endl;
            std::string tag = cfg.resolution.quantity;
            
            std::vector<double> x_value;
            std::vector<double> y_value;

            for (int i = 0; i < resolution.mean.size(); i++) {
                if (resolution.events[i] > 1) {
                    x_value.push_back(resolution.mean[i]);
                    y_value.push_back(resolution.sigma[i]);
                }
            }

            auto canvas = TCanvas(("c_res_" + tag).c_str(), ("Resolution_ " + tag).c_str(), 800, 600);

            TGraph* graph = new TGraph(x_value.size(), x_value.data(), y_value.data());

            std::string axis;
            if (tag == "pt") {
                axis = ";Mean P_t [GeV];#sigma P_t [GeV]";
            } else {
                axis = ";Mean;Sigma";
            }

            std::string title = "Resolution " + tag + axis;
            graph->SetTitle(title.c_str());
            graph->SetMarkerStyle(20);
            graph->SetMarkerSize(0.7);
            graph->SetMarkerColor(kBlue+1);
            graph->SetLineColor(kBlue+1);
            canvas.SetGrid();

            graph->Draw("AP");
            
            app->Run();
            
            delete app; 
        } else {
            std::cout << "No visualization booked." << std::endl;
        }

    }

    // ---------
    // Selection
    // ---------

    if (cfg.general.operation_mode.find("Selection") != std::string::npos) {
        try {
            ROOT::EnableImplicitMT();// MultiThread option: ON

            if (cfg.selection.dataset == "DATA") {
                dataset_tree = cfg.io.tree_data_name;
                dataset_file = cfg.io.in_data_file;
                
                std::cout << "Initializing Selection operation in DATA." << std::endl;
            } else if (cfg.selection.dataset == "MC") {
                dataset_tree = cfg.io.tree_mc_name;
                dataset_file = cfg.io.in_mc_file;
                
                std::cout << "Initializing Selection operation in MC." << std::endl;
            }

            ROOT::RDataFrame selection_data_frame(dataset_tree, dataset_file);
            
            // ----------
            // RespMatrix
            // ----------
            
            ROOT::RDF::RNode node_RM = selection_data_frame;
            
            if (cfg.selection.selection_mode.find("RespMatrix") != std::string::npos) {
                if (cfg.selection.dataset == "DATA") {
                    std::cout << "ERROR: invalid Selection dataset for Response Matrix Calculus, pls select MC dataset in Selection settup, exiting..."
                     << std::endl;
                    return 1; 
                }

                if (cfg.general.verbose) {
                    std::cout << "Executing Response Matrix selection." << std::endl;
                }
            
                node_RM = CalculateRespMatrixWrapper(node_RM, flags_RM, cuts_RM);
            }

            // -----
            // Event
            // -----

            ROOT::RDF::RNode node_event = selection_data_frame;
            
            if (cfg.selection.selection_mode.find("Event") != std::string::npos) {
                if (cfg.general.verbose) {
                    std::cout << "Executing Event selection." << std::endl;
                }

                if (cfg.selection.dataset == "DATA") {
                    node_event = ApplyValidationFilter(node_event, validation_map, "run", "luminosityBlock");
                }

                node_event = EventSelection(node_event, cfg);
            }

            // --------------
            // Output Manager
            // --------------

            OutputSelManager manager(cfg);

            if (cfg.selection.selection_mode == "RespMatrix") {
                std::cout << "Starting RespMatrix booking." << std::endl;
                manager.BookAction(node_RM, cfg);
            } else if (cfg.selection.selection_mode == "Event") {
                std::cout << "Starting Event booking." << std::endl;
                manager.BookAction(node_event, cfg);
            }

            if (cfg.general.verbose) {
                std::cout << "Starting Output Selection Manager, running..." << std::endl;
            }

            manager.Run();

            if ((cfg.selection.visual_sel) && (app != nullptr)) {                
                std::cout << "Starting visualization..." << std::endl;
                app->Run();
                
                delete app; 
            } else {
                std::cout << "No visualization booked." << std::endl;
            }

        } catch (const std::exception& except) {
            std::cout << "Error nature: " << except.what() << std::endl;
            return 1;
        }
    }

    // --------
    // Analysis
    // --------

    if (cfg.general.operation_mode == "Analysis") {
        ROOT::EnableImplicitMT();// MultiThread option: ON
        
        if (cfg.analysis.analysis_mode == "Unfold") {
            std::string tag = cfg.unfold.unfold_quantity;
            std::vector<std::unique_ptr<TCanvas>> canvas;

            // MonteCarlo DATASET
            ROOT::RDataFrame mc_frame("MC_RespMatrix_Tree", cfg.selection.o_sel_file_data);
            ROOT::RDF::RNode node_RM = mc_frame;
            
            RespMatrixHisto resp_histo = BuildRespMatrixHisto(node_RM, cfg);

            if (cfg.unfold.check_plot) {
                
                ControlHisto control_histo = BuildControlHisto(node_RM, cfg);
                
                if (visualize && app != nullptr) {
                
                    int check_control = VisualizeControlPlots(canvas, resp_histo, control_histo, tag);
                    app->Run();
                    delete app;
                    
                    return 0; 
                } else {
                    std::cout << "check_plot True but no visualization was booked, exiting..." << std::endl;
                    return 1;
                }
            }

            // DATA DATASET
            ROOT::RDataFrame data_frame("DATA_Event_Tree", cfg.selection.o_sel_file_data);
            ROOT::RDF::RNode node_event = data_frame;

            // First Event Loop on MonteCarlo
            UnfoldDensities density = CreateUnfoldDensity(resp_histo, tag);// OR HERE
            UnfoldResult result;
            
            EventSelectionHisto event_histo_struct = BuildEventSelection_Histo(node_event, cfg);

            // Signal Fitter - Second Event Loop on data

            std::unique_ptr<TH1D> event_histo;

            if (cfg.unfold.closure_test == true) {
                std::cout << "Initializing Unfold Closure test procedure." << std::endl;
                event_histo = nullptr;    
            
            } else if (cfg.unfold.bkg_subtraction == true){
                std::cout << "Initializing Unfold procedure with BKG substraction fit." << std::endl;
                event_histo = EventFit_SignalHisto_Wrapper(event_histo_struct, cfg);
            
            } else {
                std::cout << "Initializing Standard Unfold procedure." << std::endl;

                if (tag == "pt") {
                    event_histo.reset(event_histo_struct.h1_pt.GetPtr());
                } else if (tag == "y") {
                    event_histo.reset(event_histo_struct.h1_y.GetPtr());
                } else if (tag == "phis") {
                    event_histo.reset(event_histo_struct.h1_phis.GetPtr());
                }
            }
            
            // Old method
            //EventHisto event_histo = BuildEventHisto(node_event, cfg);

            if (tag == "pt") {
                result = ApplyUnfold(std::move(density.pt_unf), event_histo.get(), resp_histo.h1_pt_test.GetPtr(), 
                resp_histo.h1_pt_fake.GetPtr(), cfg, "Pt_Z0");
            
            } else if (tag == "y") {
                result = ApplyUnfold(std::move(density.y_unf), event_histo.get(), resp_histo.h1_y_test.GetPtr(), 
                resp_histo.h1_y_fake.GetPtr(), cfg, "Y_Z0");
            
            } else if (tag == "phis") {
                result = ApplyUnfold(std::move(density.phis_unf), event_histo.get(), resp_histo.h1_phis_test.GetPtr(),
                resp_histo.h1_phis_fake.GetPtr(), cfg, "Phis_Z0");
            
            } else {
                std::cout << "ERROR: invalid unfold quantity input, exiting..." << std::endl;
                return 1;
            }

            if (visualize && app != nullptr) {
                int check = VisualizeUnfoldResults(canvas, result, resp_histo, tag);
                
                if (check != 0) {
                    std::cout << "ERROR: Visualize operation failed, exiting..." << std::endl;
                    return 1;
                }
                
                app->Run();
                
                delete app;
                
                return 0;
            } else {
                std::cout << " Unfold procedure applied but no visualization was booked." << std::endl;
                return 1;
            }

        } else if (cfg.analysis.analysis_mode == "Event") {
            std::vector<std::unique_ptr<TCanvas>> canvas;

            // DATA DATASET
            ROOT::RDataFrame data_frame("DATA_Event_Tree", cfg.selection.o_sel_file_data);
            ROOT::RDF::RNode node_event = data_frame;
            EventSelectionHisto event_histo_struct = BuildEventSelection_Histo(node_event, cfg);

            std::unique_ptr<TH1D> event_histo = EventFit_SignalHisto_Wrapper(event_histo_struct, cfg);
            std::string tag = cfg.event.event_quantity;
            TH1D* histo_ev;

            if (visualize && app != nullptr) {                
                if (tag == "pt") {
                    histo_ev = event_histo_struct.h1_pt.GetPtr();
                } else if (tag == "y") {
                    histo_ev = event_histo_struct.h1_y.GetPtr();
                } else if (tag == "phis") {
                    histo_ev = event_histo_struct.h1_phis.GetPtr();
                }

                TH1D* histo = event_histo.get();
                if (histo == nullptr) {
                    std::cout << "ERROR: invalid Response Matrix histogram, exiting..." << std::endl;
                    return 1;
                }

                std::string name_c1 = "Signal yield_ " + cfg.event.event_quantity;
                auto c1 = std::make_unique<TCanvas>(name_c1.c_str(), name_c1.c_str(), 800, 600);

                std::string title_x = cfg.event.event_quantity;
                std::string title_y = cfg.event.event_quantity + " Entries";
                histo->GetXaxis()->SetTitle(title_x.c_str());
                histo->GetYaxis()->SetTitle(title_y.c_str());

                histo_ev->SetLineColor(kBlue);
                histo_ev->SetLineWidth(2);
        
                histo->SetLineColor(kRed);
                histo->SetLineWidth(2);

                double max_val = std::max(histo_ev->GetMaximum(), histo->GetMaximum());
                histo_ev->SetMaximum(max_val * 1.25);
                histo_ev->SetMinimum(0.0);

                histo_ev->Draw("E");
                histo->Draw("E SAME");
                
                c1->Update();
                canvas.push_back(std::move(c1));

                app->Run();
                
                delete app;
                
                return 0;
            } else {
                std::cout << " Unfold procedure applied but no visualization was booked." << std::endl;
                return 1;
            }
        }
    }

    return 0;
}