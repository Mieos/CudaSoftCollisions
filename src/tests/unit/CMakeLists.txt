#Unit tests

#Create the subdirectory for the binary files
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests/unit)

#Include the headers
#NONE

#Set the output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests/unit)

#Create the binary files
add_executable(uTest_extractSurfaceVolume uTest_extractSurfaceVolume.cpp)
target_link_libraries(uTest_extractSurfaceVolume cudaSoftCollisions)
