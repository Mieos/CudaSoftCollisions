#Create folder
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/generated/helpers)

#Generate const that can be used by including the pathsHelpers.hpp
set(MIEOS_HELPERS_DATAPATH ${CMAKE_SOURCE_DIR}/data )

configure_file( pathsHelpers.hpp.in ${CMAKE_BINARY_DIR}/generated/helpers/pathsHelpers.hpp )


