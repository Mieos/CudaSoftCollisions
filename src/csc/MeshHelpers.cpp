#include "csc/MeshHelpers.hpp"

#include <vtkUnstructuredGridReader.h>

//Read 3d VTK file
bool MeshHelpers::readVolumeMeshVTK(const std::string & fileName, vtkSmartPointer<vtkUnstructuredGrid> & mesh3d){

   vtkSmartPointer<vtkUnstructuredGridReader> reader3d = vtkSmartPointer<vtkUnstructuredGridReader>::New();
   reader3d->SetFileName (fileName.c_str() );
   reader3d->Update();
   mesh3d = reader3d->GetOutput();
   return true;

}

bool MeshHelpers::getSurfaceOfVolMesh(const vtkSmartPointer<vtkUnstructuredGrid> & mesh3d){

   size_t numberCells = mesh3d->GetNumberOfCells();
   
   std::cout << "DEBUG = " << numberCells << std::endl;

   //TODO
   return false;

}

