#include <ROOT/RDataFrame.hxx>
#include <iostream>
#include <string>


int main() {
    std::string tree_name = "Events";
    std::string file_name = "../data/dati0.root";

    try {
        ROOT::RDataFrame data_frame(tree_name, file_name);
        auto columns = data_frame.GetColumnNames();

        std::cout << "Number of columns finded: " << columns.size() << std::endl;

    } catch (const std::exception& except){
        std::cout << "Natura dell'errore: " << except.what() << std::endl;
        return 1;
    }
    return 0;
}