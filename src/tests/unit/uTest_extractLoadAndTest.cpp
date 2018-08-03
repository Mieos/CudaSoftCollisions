#include <iostream>
#include <string>
#include <fstream>

#include "vtkPLYWriter.h"

#include "pathsHelpers.hpp"

#include "csc/MeshStructureExtractor.hpp"
#include "csc/MeshStructureCollider.hpp"
#include "csc/MeshHelpers.hpp"

#include <map>
#include <vtkTriangle.h>
#include <vtkCellArray.h>
#include <vtkPLYWriter.h>

#include <chrono>

int main(int argc, char *argv[]){

   //Get and print the path to data
   std::string path = MIEOS_HELPERS_DATAPATH_GENERATED;

   //SavePath
   std::string savePath = MIEOS_HELPERS_DATAPATH_GENERATED;
   savePath = savePath + "/meshes/test3Dmodel.ply";

   //Icosaheron
   path = path + "/meshes/torus_vol_highres.vtk";
   std::cout << "Reading: " << path << std::endl;

   //Test model extraction
   cv::Mat matPointsTet;
   std::vector<size_t> associationResults;
   std::vector<size_t> tetSelected;
   std::vector<std::vector<size_t>> tetIdVector;
   std::cout << "Extracting model ... ";
   auto startExtraction = std::chrono::steady_clock::now(); 
   MeshStructureExtractor::extractModelFromFile(path,matPointsTet,tetIdVector,associationResults,tetSelected);
   auto endExtraction = std::chrono::steady_clock::now();
   auto diffExtraction = endExtraction - startExtraction;
   std::cout << std::chrono::duration <double, std::milli> (diffExtraction).count() << " ms" << std::endl;

   vtkSmartPointer<vtkUnstructuredGrid> mesh3d;
   MeshHelpers::readVolumeMeshVTK(path,mesh3d);
   size_t numberOfPoints = mesh3d->GetNumberOfPoints();
   cv::Mat pointsMat = cv::Mat::zeros(numberOfPoints,3,CV_32FC1);
   double interPoint[3];
   for(size_t k=0; k<numberOfPoints; k++){
      mesh3d->GetPoint(k,interPoint);
      pointsMat.at<float>(k,0)=float(interPoint[0])+1.0f;
      pointsMat.at<float>(k,1)=float(interPoint[1]);
      pointsMat.at<float>(k,2)=float(interPoint[2]);
   }

   //Create the object
   std::cout << "Creating model ... ";
   auto startCreation = std::chrono::steady_clock::now();  
   MeshStructureCollider* msc = new MeshStructureCollider(matPointsTet,tetIdVector,associationResults);
   auto endCreation = std::chrono::steady_clock::now();
   auto diffCreation = endCreation - startCreation;
   std::cout << std::chrono::duration <double, std::milli> (diffCreation).count() << " ms" << std::endl;

   //Update 
   std::cout << "Updating points ... ";
   auto startUpdate = std::chrono::steady_clock::now();
   msc->updatePointsPositions(pointsMat); 
   auto endUpdate = std::chrono::steady_clock::now();
   auto diffUpdate = endUpdate - startUpdate;
   std::cout << std::chrono::duration <double, std::milli> (diffUpdate).count() << " ms" << std::endl;

   //Collide
   std::vector<bool> collideVector;
   for(size_t k=0; k<tetSelected.size();k++){
      collideVector.push_back(false);
   }
   std::cout << "Colliding ... ";
   auto startCollide = std::chrono::steady_clock::now(); 
   msc->collide(collideVector);
   auto endCollide = std::chrono::steady_clock::now();
   auto diffCollide = endCollide - startCollide;
   std::cout << std::chrono::duration <double, std::milli> (diffCollide).count() << " ms" << std::endl;

   //Check results
   for(size_t k=0; k<collideVector.size();k++){
      if(collideVector.at(k)){
         std::cout << "Collide : " << k << std::endl;
      }
   }

   delete msc;

   std::cout << "Test : Ending.. " << std::endl << std::endl;



   return 0;
}
