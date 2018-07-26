#include <iostream>
#include <string>
#include <fstream>

#include "pathsHelpers.hpp"
#include "csc/MeshHelpers.hpp"

int main(int argc, char *argv[]){

   //Get and print the path to data
   std::string path = MIEOS_HELPERS_DATAPATH_GENERATED;
   
   //Icosaheron
   path = path + "/meshes/icosahedron_vol.vtk";
   std::cout << "Reading: " << path << std::endl;

   //Read a data file
   vtkSmartPointer<vtkUnstructuredGrid> mesh3d;
   MeshHelpers::readVolumeMeshVTK(path,mesh3d);
   std::cout << std::endl;

   //Extract surface
   MeshHelpers::getSurfaceOfVolMesh(mesh3d);
   std::cout << "Test : Ending.. " << std::endl << std::endl;

   return 0;
}
