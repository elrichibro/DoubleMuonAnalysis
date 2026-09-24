#include "Utils.h"

std::vector<double> CreateBins(int nbins, double min, double max, const std::string& distribution) {
    // Initial definition
    std::vector<double> bins;
    bins.reserve(nbins + 1);
    
    // Distribution option
    if (distribution == "linear") {
        double step = (max - min)/nbins;

        for (int i = 0; i < nbins + 1; i++){
            bins.push_back(min + i * step);
        }

    } else if (distribution == "log") {

        double log_min = std::log10(min);
        double log_max = std::log10(max);
        
        double step = (log_max - log_min)/nbins;

        for (int i = 0; i < nbins + 1; i++){
            bins.push_back(std::pow(10, log_min + i * step));

        }
    } else {
        std::cout << "ERROR: invalid distribution, exiting..." << std::endl;
        return bins;
    }

    return bins;
}