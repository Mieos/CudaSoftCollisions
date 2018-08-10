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
   path = path + "/meshes/triangularDipyramid_volume.vtk";
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
      pointsMat.at<float>(k,0)=float(interPoint[0]);
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

   //Check results (no intersection here)
   bool bugDetected= false;
   for(size_t k=0; k<collideVector.size();k++){
      if(collideVector.at(k)){
         bugDetected=true;
      }
   }

   //Print if there is a problem
   if(bugDetected){
      std::cout << "BUG detected : intersection detected while it should not" << std::endl;
   } else {
      std::cout << "No intersection detected : OK" << std::endl;
   }

   //Update point to create intersection
   float vectorChange[3];
   vectorChange[0]=0.6*(pointsMat.at<float>(2,0) - pointsMat.at<float>(0,0));
   vectorChange[1]=0.6*(pointsMat.at<float>(2,1) - pointsMat.at<float>(0,1));
   vectorChange[2]=0.6*(pointsMat.at<float>(2,2) - pointsMat.at<float>(0,2));
   pointsMat.at<float>(0,0) = pointsMat.at<float>(0,0) + vectorChange[0];
   pointsMat.at<float>(0,1) = pointsMat.at<float>(0,1) + vectorChange[1];
   pointsMat.at<float>(0,2) = pointsMat.at<float>(0,2) + vectorChange[2];

   //Update the object
   std::cout << "Updating points (creating collision)... ";
   auto startUpdate2 = std::chrono::steady_clock::now();
   msc->updatePointsPositions(pointsMat); 
   auto endUpdate2 = std::chrono::steady_clock::now();
   auto diffUpdate2 = endUpdate2 - startUpdate2;
   std::cout << std::chrono::duration <double, std::milli> (diffUpdate2).count() << " ms" << std::endl;

   //Collide
   std::cout << "Colliding ... ";
   auto startCollide2 = std::chrono::steady_clock::now(); 
   msc->collide(collideVector);
   auto endCollide2 = std::chrono::steady_clock::now();
   auto diffCollide2 = endCollide2 - startCollide2;
   std::cout << std::chrono::duration <double, std::milli> (diffCollide2).count() << " ms" << std::endl;

   //Check results (intersection here)
   bugDetected= false;
   for(size_t k=0; k<collideVector.size();k++){
      if(!collideVector.at(k)){
         bugDetected=true;
      }
   }

   //Print if there is a problem
   if(bugDetected){
      std::cout << "BUG detected : intersection not detected while it should" << std::endl;
   } else {
      std::cout << "Intersection correctly detected : OK" << std::endl;
   }

   delete msc;

   std::cout << "Test : Ending.. " << std::endl << std::endl;

   return 0;
}
