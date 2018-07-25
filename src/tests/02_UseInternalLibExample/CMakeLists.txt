#02_USEINTERNALLIBEXAMPLE

#Create the subdirectory for the binary files
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests/02_UseInternalLibExample)

#Include the headers
#NONE

#Set the output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests/02_UseInternalLibExample)

#Create the binary files
add_executable(test02Run main.cpp)

target_link_libraries(test02Run
   MIEOSPO_INTLIBRARYEXAMPLE_LIB
   )
