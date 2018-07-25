#Minimum cmake version asked
cmake_minimum_required(VERSION 2.8)

#Project Name
project(Mieos-Project-Organizer-Example)

#Use C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#Other Flags you want to add to the project
#For instance here, allow pthread
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

#Options
#option(myOption "myOption" ON)

##########
#Packages#
##########

######
#Srcs#
######

##
#Pkg config allow you to find unconventionnal packages : uncomment the next line if you need it
##
#find_package(PkgConfig REQUIRED)

##
#Other packages you want to include with the classic method
##
#find_package(myPackage REQUIRED)

##
#Other packages you want to include with PkgConfig
##
#pkg_search_module(MyPackage REQUIRED myPackagev42)

##
#External packages
##
add_subdirectory(external)

##########
#Includes#
##########

##
#Classic packages
##
#if( MYPACKAGE_FOUND )
#   include_directories(${MYPACKAGE_INCLUDE_DIRS})
#endif()

##
#External packages
##
include_directories(
   external/externalLibraryExample/include
   )

#LIBS
set(ALL_LIBS
   MIEOSPO_EXTLIBEXAMPLE_LIB
   )

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
#internalLibraryExample
add_subdirectory(src/internalLibraryExample)
include_directories(include/internalLibraryExample)
#Tests
add_subdirectory(src/tests)
