#Minimum cmake version asked
cmake_minimum_required(VERSION 2.8)

#Project Name
project(CudaSoftCollisions)

#Use C++11
set(CMAKE_CXX_STANDARD 11)

##########
#Packages#
##########

find_package(OpenCV REQUIRED)
find_package(CUDA REQUIRED)
find_package(VTK REQUIRED)

######
#Srcs#
######

##########
#Includes#
##########
#Does not include vtk headers here because vtk does not like them
include_directories(include)

###################
#Generated helpers#
###################
#Generated helpers such as paths and others stuffs that can be used in cpp afterward
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
add_subdirectory(src/csc/cuda_src) #Cuda does not like vtk headers..
include(${VTK_USE_FILE})
add_subdirectory(src/csc)

#Tests
add_subdirectory(src/tests)
