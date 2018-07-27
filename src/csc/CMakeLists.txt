#Lib

set(sudaSoftCollisions_src
   MeshStructureExtractor.cpp
   MeshHelpers.cpp
   )

add_library(cudaSoftCollisions SHARED ${sudaSoftCollisions_src})
target_link_libraries(cudaSoftCollisions ${VTK_LIBRARIES}  ${OpenCV_LIBS})
export(PACKAGE cudaSoftCollisions)



