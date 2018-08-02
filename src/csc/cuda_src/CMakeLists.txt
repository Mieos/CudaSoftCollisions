#Lib

set(cudaSoftCollisions_src_cuda
   MeshStructureCollider.cu
   )

cuda_add_library(cudaSoftCollisions_cuda ${cudaSoftCollisions_src_cuda})
target_link_libraries(cudaSoftCollisions_cuda ${OpenCV_LIBS})

export(PACKAGE cudaSoftCollisions_cuda)



