#INTERNALLIBRARYEXAMPLE

##Include the headers
include_directories(${CMAKE_SOURCE_DIR}/include/internalLibraryExample)

#Set all the sources
set(MIEOSPO_INTLIBRARYEXAMPLE_SOURCES_LIB
   MieosProjectOrganiserIntLib.cpp
   )

#We create the library
add_library(MIEOSPO_INTLIBRARYEXAMPLE_LIB ${MIEOSPO_INTLIBRARYEXAMPLE_SOURCES_LIB} )

export(PACKAGE MIEOSPO_INTLIBRARYEXAMPLE_LIB)
