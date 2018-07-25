#Lib

set(sudaSoftCollisions_src
   MeshStructureExtractor.cpp
   )

add_library(sudaSoftCollisions SHARED ${sudaSoftCollisions_src})
export(PACKAGE sudaSoftCollisions)



