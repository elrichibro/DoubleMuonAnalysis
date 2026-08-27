cmake_minimum_required(VERSION 3.10)

project(DoubleMuonAnalysis CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ROOT REQUIRED COMPONENTS ROOTDataFrame)
find_package(nlohmann_json REQUIRED)

include_directories(include)

set(COMMON_SRCS
    src/Config.cpp
    src/Filters.cpp
)

set(COMMON_LIBS
    ROOT::ROOTDataFrame
    nlohmann_json::nlohmann_json
)

add_executable(analyze src/main.cpp ${COMMON_SRCS})
target_link_libraries(analyze PRIVATE ${COMMON_LIBS})

add_executable(analyze_mc src/mainMC.cpp ${COMMON_SRCS})
target_link_libraries(analyze_mc PRIVATE ${COMMON_LIBS})

add_executable(analyze_merge src/main_merge.cpp ${COMMON_SRCS})
target_link_libraries(analyze_merge PRIVATE ${COMMON_LIBS})