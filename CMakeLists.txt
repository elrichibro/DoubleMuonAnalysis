cmake_minimum_required(VERSION 3.10)

project(DoubleMuonAnalysis CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ROOT REQUIRED COMPONENTS ROOTDataFrame)
find_package(nlohmann_json REQUIRED)

include_directories(include)

add_executable(analyze src/main.cpp src/Config.cpp)

target_link_libraries(analyze PRIVATE ROOT::ROOTDataFrame nlohmann_json::nlohmann_json)
