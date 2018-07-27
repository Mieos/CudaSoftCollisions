#include "csc/MeshStructureExtractor.hpp"
#include "csc/MeshHelpers.hpp"

bool MeshStructureExtractor::extractModelFromFile(std::string fileName){

   vtkSmartPointer<vtkUnstructuredGrid> msh;
   if(!MeshHelpers::readVolumeMeshVTK(fileName,msh)){
      return false;
   }

   return extractModelFromMesh(msh);

}

bool MeshStructureExtractor::extractModelFromMesh(const vtkSmartPointer<vtkUnstructuredGrid> & mesh3d){

   //First get the surface of the model
   vtkSmartPointer<vtkPolyData> surfPolyData;
   MeshHelpers::getSurfaceOfVolMesh(mesh3d, surfPolyData);

   return false; //TODO

}