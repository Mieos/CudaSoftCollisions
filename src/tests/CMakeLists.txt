#TESTS

#Create folder
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests)

#Add subdirectories
add_subdirectory(01_UseExternalLibExample)
add_subdirectory(02_UseInternalLibExample)
add_subdirectory(03_UseHelpers)
