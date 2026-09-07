cmake_minimum_required(VERSION 3.10)

project(DoubleMuonAnalysis CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ROOT REQUIRED COMPONENTS ROOTDataFrame)
find_package(nlohmann_json REQUIRED)
find_package(ROOT REQUIRED COMPONENTS Core Hist Tree RIO MathCore RooFit RooFitCore Minuit)

include_directories(include)

set(COMMON_SRCS
    src/Config.cpp
    src/Filters.cpp
    src/Checks.cpp
    src/Manager.cpp
    src/AnalysisTools.cpp
)

set(COMMON_LIBS
    ROOT::ROOTDataFrame
    nlohmann_json::nlohmann_json
    ROOT::Core 
    ROOT::Hist 
    ROOT::Tree 
    ROOT::RIO 
    ROOT::MathCore 
    ROOT::RooFit 
    ROOT::RooFitCore
    ROOT::Minuit
)

#add_executable(analyze src/main.cpp ${COMMON_SRCS})
#target_link_libraries(analyze PRIVATE ${COMMON_LIBS})

add_executable(analyze_mc src/mainMC.cpp ${COMMON_SRCS})
target_link_libraries(analyze_mc PRIVATE ${COMMON_LIBS})