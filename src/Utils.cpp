#include "Utils.h"

std::vector<float> CreateBins(int nbins, float min, float max, const std::string& distribution, float split) {
    std::vector<float> bins;
    bins.reserve(nbins + 1);
    
    if (distribution == "uniform") {
        float step = (max - min)/nbins;

        for (int i = 0; i <= nbins; i++){
            bins.push_back(min + i * step);
        }

    } else if (distribution == "log") {

        float log_min = std::log10(min);
        float log_max = std::log10(max);
        
        float step = (log_max - log_min)/nbins;

        for (int i = 0; i <= nbins; i++){
            bins.push_back(std::pow(10, log_min + i * step));

        }
    } else {
        std::cout << "ERROR: invalid distribution, exiting..." << std::endl;
        return bins;
    }

    return bins;
}