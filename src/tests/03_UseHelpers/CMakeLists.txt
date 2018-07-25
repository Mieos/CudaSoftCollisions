#03_USEHELPERS

#Create the subdirectory for the binary files
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests/03_UseHelpers)

#Include the headers
#NONE

#Set the output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests/03_UseHelpers)

#Create the binary files
add_executable(test03Run main.cpp)
