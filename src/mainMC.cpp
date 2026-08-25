#include <ROOT/RDataFrame.hxx>

#include <iostream>
#include <string>

#include "Config.h"
#include "Filters.h"

int main(int argc, char* argv[]) {

    // Variables needed
    config_struct cfg;
    std::string json_path = "";
    int verbose = 0;
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg.rfind(".json") != std::string::npos) {
            json_path = arg;
        } else if (arg == "--verbose" || arg == "-v") {
            verbose = 1;
        } else {
            std::cout << "ERROR: invalid input command, please try again, exinting.\n" << std::endl;
            return 1;
        }
    }
/*
    if (verbose) {
        Verbose_config(cfg);
    }
*/
    if (Configure(cfg, json_path) != 0) {
        std::cout << "Configurations fails, check needed, exiting." << std::endl;
        return 1;
    }

    if (verbose) {
        Verbose_config(cfg);
    }
    
    try {
        // Trying to Enable Implicit MT
        ROOT::EnableImplicitMT();

        if (verbose){ 
            std::cout << "RDataFrame object created, starting analysis ..." << std::endl;
        }

        ROOT::RDataFrame data_frame(cfg.io.tree_mc_name, cfg.io.in_mc_file);

        // Validation run filter
        auto node_is_Z0 = data_frame
        .Define("is_good_Z0", is_Good_Z0, {"GenPart_pdgId", "GenPart_statusFlags"})
        .Filter("Sum(is_good_Z0) == 1", "1. True Z0");

        auto cut_report = node_is_Z0.Report();

        cut_report->Print();

    } catch (const std::exception& except) {
        std::cerr << "Error nature: " << except.what() << std::endl;
        return 1;
    }

    return 0;
}