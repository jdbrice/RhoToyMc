cmake_minimum_required(VERSION 3.10)

project(toy)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(toy Engine.cpp)

target_include_directories(toy PUBLIC "${PROJECT_SOURCE_DIR}/" "${PROJECT_SOURCE_DIR}/XmlConfig" "${PROJECT_SOURCE_DIR}/RootAna" "${PROJECT_SOURCE_DIR}/TaskEngine" $ENV{ROOTSYS}/include )

add_subdirectory(XmlConfig)
include_directories(XmlConfig)
link_directories(XmlConfig)
add_subdirectory(TaskEngine)
include_directories(TaskEngine)
link_directories(TaskEngine)
add_subdirectory(RooPlotLib)
include_directories(RooPlotLib)
link_directories(RooPlotLib)

add_subdirectory(RootAna)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)



target_link_libraries(toy PUBLIC XmlConfig TaskEngine RootAna RooPlotLib ${ROOT_LIBRARIES})



execute_process(COMMAND root-config --cflags OUTPUT_VARIABLE ROOTLINKFLAGS)

# target_link_options(toy PUBLIC "$ROOTLINKFLAGS")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}" "${ROOTLINKFLAGS}")