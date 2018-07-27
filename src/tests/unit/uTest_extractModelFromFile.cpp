#include <iostream>
#include <string>
#include <fstream>

#include "vtkPLYWriter.h"

#include "pathsHelpers.hpp"

#include "csc/MeshStructureExtractor.hpp"

int main(int argc, char *argv[]){

   //Get and print the path to data
   std::string path = MIEOS_HELPERS_DATAPATH_GENERATED;
   
   //Icosaheron
   path = path + "/meshes/torus_vol.vtk";
   std::cout << "Reading: " << path << std::endl;

   //Test model extraction
   MeshStructureExtractor::extractModelFromFile(path);
   std::cout << "Test : Ending.. " << std::endl << std::endl;


   
   return 0;
}
