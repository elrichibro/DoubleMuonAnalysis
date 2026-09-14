cmake_minimum_required(VERSION 3.14)

project(DoubleMuonAnalysis CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ROOT REQUIRED COMPONENTS ROOTDataFrame Core Hist Tree RIO MathCore RooFit RooFitCore Minuit)
find_package(nlohmann_json REQUIRED)

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

include(FetchContent)

FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(analyze_mc src/mainMC.cpp ${COMMON_SRCS})
target_link_libraries(analyze_mc PRIVATE ${COMMON_LIBS})

# Eseguibile per la suite di Unit Test
add_executable(run_tests 
    tests/TestUtils.cpp
    tests/TestChecks.cpp
    ${COMMON_SRCS}
)

target_link_libraries(run_tests PRIVATE gtest_main ${COMMON_LIBS})

include(GoogleTest)
gtest_discover_tests(run_tests)