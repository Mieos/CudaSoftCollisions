#TESTS

#Create folder
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests)

#Add subdirectories
add_subdirectory(03_UseHelpers)
add_subdirectory(unit)
