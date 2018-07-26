#Minimum cmake version asked
cmake_minimum_required(VERSION 2.8)

#Project Name
project(CudaSoftCollisions)

#Use C++11
set(CMAKE_CXX_STANDARD 11)

#Other Flags you want to add to the project
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

#Options
#option(myOption "myOption" ON)

##########
#Packages#
##########

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

######
#Srcs#
######

##
#External packages
##
add_subdirectory(external)

##########
#Includes#
##########

include_directories(include)

###################
#Generated helpers#
###################
#Generated helpers such as paths and others stuffs that can be used in cpp afterward
#Usefull when you need an absolute path
add_subdirectory(helpers)
include_directories(${CMAKE_BINARY_DIR}/generated/helpers )

#########
#Outputs#
#########
#Set outputs dir
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

################
#Subdirectories#
################

#Lib
add_subdirectory(src/csc)

#Tests
add_subdirectory(src/tests)
