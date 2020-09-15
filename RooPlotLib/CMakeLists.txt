cmake_minimum_required(VERSION 3.10)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)

add_definitions( "-DJDB_LOG_LEVEL=0")

add_library(RooPlotLib SHARED RooPlotLib.cpp )
target_link_libraries(RooPlotLib PUBLIC ${ROOT_LIBRARIES} XmlConfig)
target_include_directories(RooPlotLib PUBLIC "$ENV{ROOTSYS}/include" )
