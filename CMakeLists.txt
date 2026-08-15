cmake_minimum_required(VERSION 3.10)
project(DoubleMuonAnalysis CXX)
set(CMAKE_CXX_STANDARD 17)

find_package(ROOT REQUIRED COMPONENTS ROOTDataFrame)

add_executable(analyze src/main.cpp)

target_link_libraries(analyze PRIVATE ROOT::ROOTDataFrame)
