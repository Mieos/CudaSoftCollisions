#include <iostream>
#include <string>
#include <fstream>

#include "vtkPLYWriter.h"

#include "pathsHelpers.hpp"
#include "csc/MeshHelpers.hpp"

int main(int argc, char *argv[]){

   //Get and print the path to data
   std::string path = MIEOS_HELPERS_DATAPATH_GENERATED;
   
   //Icosaheron
   path = path + "/meshes/torus_vol.vtk";
   std::cout << "Reading: " << path << std::endl;

   //Read a data file
   vtkSmartPointer<vtkUnstructuredGrid> mesh3d;
   MeshHelpers::readVolumeMeshVTK(path,mesh3d);
   std::cout << std::endl;

   //Extract surface
   vtkSmartPointer<vtkPolyData> surfPolyData;
   MeshHelpers::getSurfaceOfVolMesh(mesh3d, surfPolyData);
   std::cout << "Test : Ending.. " << std::endl << std::endl;

   //Save result (change that later TODO)
   std::string pathSave = MIEOS_HELPERS_DATAPATH_GENERATED;
   pathSave = pathSave + "/meshes/test.ply";
   vtkSmartPointer<vtkPLYWriter> plyWriter = vtkSmartPointer<vtkPLYWriter>::New();
   plyWriter->SetFileName(pathSave.c_str());
   plyWriter->SetFileTypeToASCII();

#if VTK_MAJOR_VERSION <= 5
   plyWriter->SetInput(surfPolyData);
#else
   plyWriter->SetInputData(surfPolyData);
#endif
   plyWriter->Write();
   
   return 0;
}
