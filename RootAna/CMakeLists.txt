cmake_minimum_required(VERSION 3.10)

message( STATUS "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS} )
#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)

add_library(RootAna SHARED HistoAnalyzer.cpp HistoBook.cpp Reporter.cpp )
target_link_libraries(RootAna PUBLIC ${ROOT_LIBRARIES} XmlConfig RooPlotLib)

add_definitions( "-DJDB_LOG_LEVEL=0")

target_include_directories(RootAna PUBLIC "$ENV{ROOTSYS}/include" )
