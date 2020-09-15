cmake_minimum_required(VERSION 3.10)
add_library(TaskEngine SHARED TaskEngine.cpp TaskFactory.cpp )
target_link_libraries(TaskEngine PUBLIC XmlConfig)
add_definitions( "-DJDB_LOG_LEVEL=0")
