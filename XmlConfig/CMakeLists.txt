cmake_minimum_required(VERSION 3.10)

# project(XMLConfig)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_library(XmlConfig SHARED XmlConfig.cpp Utils.cpp XmlString.cpp Logger.cpp )

# set_property(TARGET XMLConfig PROPERTY CXX_STANDARD 11)


add_definitions( "-DJDB_LOG_LEVEL=0")

file( GLOB all_headers "*.h")
install(TARGETS XmlConfig DESTINATION lib)
install(FILES ${all_headers} DESTINATION include/XmlConfig)