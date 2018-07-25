#01_USEEXTERNALLIBEXAMPLE

#Create the subdirectory for the binary files
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests/01_UseExternalLibExample)

#Include the headers
#NONE

#Set the output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests/01_UseExternalLibExample)

#Create the binary files
add_executable(test01Run main.cpp)

target_link_libraries(test01Run
   MIEOSPO_EXTLIBEXAMPLE_LIB
   )
