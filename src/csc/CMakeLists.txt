#Lib

set(cudaSoftCollisions_src
   MeshStructureExtractor.cpp
   MeshHelpers.cpp
   )

add_library(cudaSoftCollisions SHARED ${cudaSoftCollisions_src})
target_link_libraries(cudaSoftCollisions ${VTK_LIBRARIES}  ${OpenCV_LIBS} cudaSoftCollisions_cuda)
export(PACKAGE cudaSoftCollisions)



