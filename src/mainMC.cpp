#include <ROOT/RDataFrame.hxx>

#include <iostream>
#include <string>

#include <TApplication.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMarker.h>
#include "TEfficiency.h"

#include "Config.h"
#include "Filters.h"
#include "Utils.h"
#include "Checks.h"
#include "Manager.h"
#include "AnalysisTools.h"
#include "Unfold.h"

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
    const flags_config flags_TP = cfg.flag_TP;
    const cuts_config cuts_TP = cfg.cut_TP;

    const flags_config flags_RM = cfg.flag_RM;
    const cuts_config cuts_RM = cfg.cut_RM;

    validation_type validation_map = Validation_load(cfg.io.val_file);
    std::cout << "Validation Map created." << std::endl;

    std::string dataset_tree = "";
    std::string dataset_file = "";
    
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

            // -------------
            // Tag and Probe
            // -------------
            
            ROOT::RDF::RNode node_TP = selection_data_frame;
            
            if (cfg.selection.selection_mode.find("TagAndProbe") != std::string::npos) {
                
                if (cfg.general.verbose) {
                    std::cout << "Executing TagAndProbe selection." << std::endl;
                }

                node_TP = CalculateTagAndProbeWrapper(node_TP, validation_map, cfg.selection, flags_TP, cuts_TP);
            }

            // -----
            // Event
            // -----

            ROOT::RDF::RNode node_event = selection_data_frame;
            
            if (cfg.selection.selection_mode.find("Event") != std::string::npos) {
                if (cfg.selection.dataset == "MC") {
                    std::cout << "ERROR: invalid Selection dataset for Event selection, pls select DATA dataset in this Selection settup, exiting..."
                     << std::endl;
                    return 1; 
                }

                if (cfg.general.verbose) {
                    std::cout << "Executing Event selection." << std::endl;
                }
                node_event = ApplyValidationFilter(node_event, validation_map, "run", "luminosityBlock");
                
                node_event = EventSelection(node_event, cfg);
            }

            // --------------
            // Output Manager
            // --------------

            OutputSelManager manager(cfg);

            if (cfg.selection.selection_mode == "TagAndProbe") {
                std::cout << "Starting TagAndProbe booking." << std::endl;
                manager.BookAnalysis(node_TP, cfg);
            } else if (cfg.selection.selection_mode == "RespMatrix") {
                std::cout << "Starting RespMatrix booking." << std::endl;
                manager.BookAnalysis(node_RM, cfg);
            } else if (cfg.selection.selection_mode == "Event") {
                std::cout << "Starting Event booking." << std::endl;
                manager.BookAnalysis(node_event, cfg);
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
    // Template
    // --------

    if (cfg.general.operation_mode.find("Template") != std::string::npos) {
        ROOT::EnableImplicitMT();
        
        if (verbose) {
            std::cout << "Initializating Template Operation Mode..." << std::endl;
        }

        if (cfg.templ.template_type.find("UNBINNED") != std::string::npos) {
            if (cfg.templ.dataset == "DATA") {
                if (cfg.general.verbose) {
                    std::cout << "Transforming Selected: TagAndProbe data into -> UNBINNED data from -> DATA dataset" << std::endl;
                }
                std::string tree = cfg.templ.dataset + "_TagAndProbe_Tree";
                
                ROOT::RDataFrame data_frame(tree, cfg.selection.o_sel_file_data);
                ROOT::RDF::RNode node = data_frame;
            
                int chec_roll = UnbinnedTemplateMaker(node, cfg);  
            }

            if (cfg.templ.dataset == "MC") {
                if (cfg.general.verbose) {
                    std::cout << "Transforming Selected: TagAndProbe data into -> UNBINNED data from -> MC dataset" << std::endl;
                }
                std::string tree = cfg.templ.dataset + "_TagAndProbe_Tree";
                
                ROOT::RDataFrame data_frame(tree, cfg.selection.o_sel_file_data);
                ROOT::RDF::RNode node = data_frame;
            
                int chec_roll = UnbinnedTemplateMaker(node, cfg);  
            }
        }
        
        if (cfg.templ.template_type.find("BINNED") != std::string::npos) {
            if ((cfg.templ.dataset == "DATA")) {
                if (cfg.general.verbose) {
                    std::cout << "Transforming Selected: TagAndProbe data into -> BINNED data from -> DATA dataset" << std::endl;
                }
                std::string tree = "DATA_TagAndProbe_Tree";
                
                ROOT::RDataFrame data_frame(tree, cfg.selection.o_sel_file_data);
                ROOT::RDF::RNode node = data_frame;
                
                int chec_maker = BinnedTemplateMaker(node, cfg, 1);
            }

            if ((cfg.templ.dataset == "MC")) {
                if (cfg.general.verbose) {
                    std::cout << "Transforming Selected: TagAndProbe data into -> BINNED data from -> MC dataset" << std::endl;
                }
                std::string tree = "MC_TagAndProbe_Tree";
                
                ROOT::RDataFrame data_frame(tree, cfg.selection.o_sel_file_data);
                ROOT::RDF::RNode node = data_frame;
                
                int check_maker = BinnedTemplateMaker(node, cfg, 2);
            }
        }

        if (cfg.general.verbose) {
            std::cout << "Template successfully writted into: " << cfg.templ.o_template_file_data 
            << " -> [ " << cfg.templ.bins_settup << " ]" <<  std::endl;
        }
    }

    // --------
    // Analysis
    // --------

    if (cfg.general.operation_mode == "Analysis") {
        if (cfg.analysis.analysis_mode == "TagAndProbe_DATA") {
            
            try {     
                ROOT::EnableImplicitMT();

                TFile o_template_file(cfg.templ.o_template_file_data.c_str(), "UPDATE");// Template file -> read intput

                std::vector<Template_RooF> template_container;
                std::vector<FitResult> fit_results;

                if (cfg.general.verbose) {
                    std::cout << "Loading template..." << std::endl;
                }

                if (LoadBinnedTemplate(cfg, template_container) != 0) {
                    std::cout << "ERROR: Load operations fails, exiting." << std::endl;
                    return 1;
                }
                
                TFile o_fit_file(cfg.analysis.o_fit_file.c_str(), "UPDATE");// Current writing file
                o_fit_file.cd();

                if (cfg.general.verbose) {
                    std::cout << "Starting Fit operation..." << std::endl;
                }

                int check_fit = EfficiencyFitter(template_container, cfg, fit_results, &o_fit_file);

                if (check_fit != 0) {
                    std::cout << "ERROR: Fit operation fails." << std::endl;
                    return 1;
                }
                
                if (cfg.general.verbose) {
                    std::cout << "Fit operation finished." << std::endl;
                }

                std::vector<std::string> booked_values = {"efficiency", "n_tot", "fit_status", "mu", "sigma", "lambda_pass", "lambda_fail"};
                int check = SaveMapFittedValues(&o_fit_file, fit_results, cfg, booked_values);

                if (check != 0) {
                    std::cout << "ERROR: Save operation fails, exiting..." << std::endl;
                    return 0;
                }
                
                int i = 1;
                
                if (cfg.general.verbose) {
                    std::cout << "Fit results: " << std::endl;
                    for (const auto& it : fit_results) {
                        std::cout << "" << std::endl;
                        std::cout << "Fit number: " << i << ", status: " << it.fit_status << std::endl; 
                        std::cout << "" << std::endl;
                        
                        std::cout << "    Efficiency: " << it.efficiency << " +- " << it.efficiency_err << std::endl;
                        std::cout << "    Total signal events: " << it.n_tot << " +- " << it.n_tot_err << std::endl;

                        std::cout << "    Mean: " << it.mu << " +- " << it.mu_err << std::endl;
                        std::cout << "    Sigma: " << it.sigma << " +- " << it.sigma_err << std::endl;
                        
                        std::cout << "    Lambda pass: " << it.lambda_pass << " +- " << it.lambda_pass_err << std::endl;
                        std::cout << "    Lambda fail: " << it.lambda_fail << " +- " << it.lambda_fail_err << std::endl;

                        i++;
                    }
                }
                
            } catch (const std::exception& except) {
                std::cerr << "Error nature: " << except.what() << std::endl;
                return 1;
            }
            
        } else if (cfg.analysis.analysis_mode == "TagAndProbe_MC") {
            
            TFile o_fit_file(cfg.analysis.o_fit_file.c_str(), "UPDATE");// Current writing file

            o_fit_file.cd();
            std::string tree = "MC_TagAndProbe_Tree";
            
            ROOT::RDataFrame data_frame(tree, cfg.selection.o_sel_file_data);
            
            if (verbose) { 
                std::cout << "RDataFrame object created, unpacking " << "MC_" + cfg.selection.selection_mode + "_Tree"
                << " from " <<  cfg.selection.o_sel_file_data << " file, starting analysis ..." << std::endl;
            }

            OutputSelManager Analysis_manager(cfg);

            Analysis_manager.BookAnalysis(data_frame, cfg);
            
            Analysis_manager.Run();
            
            if (visualize && app != nullptr) {
                std::cout << "Initializing visualization ..." << std::endl;
                app->Run();
                
                delete app; 
            } else {
                std::cout << "No visualization booked.\n" << std::endl;
            }
        
        } else if (cfg.analysis.analysis_mode == "Unfold") {
            
            // MonteCarlo DATASET
            ROOT::RDataFrame mc_frame("MC_RespMatrix_Tree", cfg.selection.o_sel_file_data);
            ROOT::RDF::RNode node_RM = mc_frame;
            
            // DATA DATASET
            ROOT::RDataFrame data_frame("DATA_Event_Tree", cfg.selection.o_sel_file_data);
            ROOT::RDF::RNode node_event = data_frame;

            // ----------
            // MonteCarlo
            // ----------

            // Insert FIT DATA procedure HERE !!!

            // First Event Loop on MonteCarlo
            RespMatrixHisto resp_histo = BuildRespMatrixHisto(node_RM, cfg);
            ControlHisto control_histo = BuildControlHisto(node_RM, cfg);// HERE

            UnfoldDensities density = CreateUnfoldDensity(resp_histo);// OR HERE

            // ----
            // DATA
            // ----

            EventHisto event_histo = BuildEventHisto(node_event, cfg);

            // Second Event Loop on DATA
            UnfoldResult result;

            if (cfg.unfold.check_plot == true) {
                std::cout << "Starting Check." << std::endl;
            
            } else if (cfg.unfold.unfold_quantity == "pt") {
                result = ApplyUnfold(std::move(density.pt_unf), event_histo.h1_pt.GetPtr(), resp_histo.h1_pt_test.GetPtr(), 
                resp_histo.h1_pt_fake.GetPtr(), cfg, "Pt_Z0");
            
            } else if (cfg.unfold.unfold_quantity == "y") {
                result = ApplyUnfold(std::move(density.y_unf), event_histo.h1_y.GetPtr(), resp_histo.h1_y_test.GetPtr(),
                resp_histo.h1_y_fake.GetPtr(), cfg, "Y_Z0");
            
            } else if (cfg.unfold.unfold_quantity == "phis") {
                result = ApplyUnfold(std::move(density.phis_unf), event_histo.h1_phis.GetPtr(), resp_histo.h1_phis_test.GetPtr(),
                resp_histo.h1_phis_fake.GetPtr(), cfg, "Phis_Z0");
            
            } else {
                std::cout << "ERROR: invalid unfold quantity input, exiting..." << std::endl;
                return 1;
            }


            std::vector<std::unique_ptr<TCanvas>> canvas;

            if (visualize && app != nullptr) {
                
                int check_control = VisualizeControlPlots(canvas, resp_histo, control_histo, cfg.unfold.unfold_quantity);

                if (cfg.unfold.check_plot == false) {
                    int check = VisualizeUnfoldResults(canvas, result, resp_histo,cfg.unfold.unfold_quantity);
                    if (check != 0) {
                        std::cout << "ERROR: Visualize operation failed, exiting..." << std::endl;
                        return 1;
                    }
                }
                app->Run();
                
                delete app; 
            } else {
                std::cout << "No visualization booked.\n" << std::endl;
            }

        } else if (cfg.analysis.analysis_mode == "Event") {
            // DATA DATASET
            ROOT::RDataFrame data_frame("DATA_Event_Tree", cfg.selection.o_sel_file_data);
            ROOT::RDF::RNode node_event = data_frame;

            EventHisto event_histo = BuildEventHisto(node_event, cfg);

            std::vector<std::unique_ptr<TH1D>> event_container = PrepareEventFit(event_histo, "pt");

            std::string file_name = "../output/Event_Fit_" + cfg.unfold.unfold_quantity + ".root"; 
            TFile o_fit_file(file_name.c_str(), "UPDATE");
            if (o_fit_file.IsZombie()) {
                std::cout << "ERROR: invalid output file, exiting..." << file_name << std::endl;
                return 1;
            }
            o_fit_file.cd();

            std::string dir_name = "fits_" + cfg.unfold.unfold_quantity;
            TDirectory* fit_dir = o_fit_file.GetDirectory(dir_name.c_str());

            if (!fit_dir) {
               fit_dir = o_fit_file.mkdir(dir_name.c_str());
            }

            std::vector<EventFitResult> event_results = EventFitWrapper(event_container, fit_dir, , const std::string& tag);

            o_fit_file.cd();
            o_fit_file.Write();
            o_fit_file.Close();
        }
    }

    return 0;
}