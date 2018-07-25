#Minimum cmake version asked
cmake_minimum_required(VERSION 2.8)

#Project Name
project(Mieos-Project-Organizer-Example-External-Library-Example)

#Use C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

##Include the headers
include_directories(include)

#Set all the sources
set(MIEOSPO_EXTLIBEXAMPLE_SOURCES_LIB
   src/MieosProjectOrganiserExtLib.cpp 
)

add_library(MIEOSPO_EXTLIBEXAMPLE_LIB ${MIEOSPO_EXTLIBEXAMPLE_SOURCES_LIB})

export(PACKAGE MIEOSPO_EXTLIBEXAMPLE_LIB)
