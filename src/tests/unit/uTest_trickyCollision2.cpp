#include <iostream>
#include <string>
#include <fstream>

#include "pathsHelpers.hpp"

#include "csc/MeshStructureExtractor.hpp"
#include "csc/MeshStructureCollider.hpp"
#include "csc/MeshHelpers.hpp"

#include <vtkPLYReader.h>
#include <chrono>

int main(int argc, char *argv[]){

   //Get and print the path to data
   std::string path = MIEOS_HELPERS_DATAPATH_GENERATED;

   //Not intersecting
   path = path + "/meshes/trickyTets2.ply";

   vtkSmartPointer<vtkPLYReader> reader = vtkSmartPointer<vtkPLYReader>::New(); 
   reader->SetFileName ( path.c_str() );
   reader->Update();
   vtkSmartPointer<vtkPolyData> mesh3D = reader->GetOutput();

   //Copy points
   size_t numPoints = mesh3D->GetNumberOfPoints();
  
   cv::Mat cvPoints = cv::Mat(numPoints,3,CV_32FC1);
   double interPoint[3];
   //Not intersecting
   for(size_t k=0; k<mesh3D->GetNumberOfPoints(); k++){
      mesh3D->GetPoint(k,interPoint);
      cvPoints.at<float>(k,0)=float(interPoint[0]);
      cvPoints.at<float>(k,1)=float(interPoint[1]);
      cvPoints.at<float>(k,2)=float(interPoint[2]);
   }

   std::vector<std::vector<size_t>> tetIdVector; 
   float v1[3], v2[3], v3[3], v1Vv2[3];
   //TET1
   for(size_t k=0; k<3; k++){
      v1[k] = cvPoints.at<float>(1,k) - cvPoints.at<float>(0,k);
      v2[k] = cvPoints.at<float>(2,k) - cvPoints.at<float>(0,k);
      v3[k] = cvPoints.at<float>(3,k) - cvPoints.at<float>(0,k);
   }
   v1Vv2[0] = v1[1]*v2[2] - v1[2]*v2[1];
   v1Vv2[1] = v1[2]*v2[0] - v1[0]*v2[2];
   v1Vv2[2] = v1[0]*v2[1] - v1[1]*v2[0];
   float dotT = v1Vv2[0]*v3[0] + v1Vv2[1]*v3[1] + v1Vv2[2]*v3[2];
   if(dotT<0){
      std::vector<size_t> interVect;
      interVect.push_back(0);
      interVect.push_back(1);
      interVect.push_back(2);
      interVect.push_back(3);
      tetIdVector.push_back(interVect);
   } else {
      std::vector<size_t> interVect;
      interVect.push_back(0);
      interVect.push_back(2);
      interVect.push_back(1);
      interVect.push_back(3);
      tetIdVector.push_back(interVect);
   }
   //TET2
   for(size_t k=0; k<3; k++){
      v1[k] = cvPoints.at<float>(5,k) - cvPoints.at<float>(4,k);
      v2[k] = cvPoints.at<float>(6,k) - cvPoints.at<float>(4,k);
      v3[k] = cvPoints.at<float>(7,k) - cvPoints.at<float>(4,k);
   }
   v1Vv2[0] = v1[1]*v2[2] - v1[2]*v2[1];
   v1Vv2[1] = v1[2]*v2[0] - v1[0]*v2[2];
   v1Vv2[2] = v1[0]*v2[1] - v1[1]*v2[0];
   dotT = v1Vv2[0]*v3[0] + v1Vv2[1]*v3[1] + v1Vv2[2]*v3[2];
   if(dotT<0){
      std::vector<size_t> interVect;
      interVect.push_back(4);
      interVect.push_back(5);
      interVect.push_back(6);
      interVect.push_back(7);
      tetIdVector.push_back(interVect);
   } else {
      std::vector<size_t> interVect;
      interVect.push_back(4);
      interVect.push_back(6);
      interVect.push_back(5);
      interVect.push_back(7);
      tetIdVector.push_back(interVect);
   }


   //Vector association (identity)
   std::vector<size_t> associationResults;
   for(size_t k=0; k<numPoints ; k++){
      associationResults.push_back(k);
   }

   for(size_t k=0; k<cvPoints.rows; k++){
      //std::cout << "PX = " << cvPoints.at<float>(k,0) << " , PY = " << cvPoints.at<float>(k,1) << " , PZ = " << cvPoints.at<float>(k,2) << std::endl;
   }

   //Create the object
   std::cout << "Creating model ... ";
   auto startCreation = std::chrono::steady_clock::now();  
   MeshStructureCollider* msc = new MeshStructureCollider(cvPoints,tetIdVector,associationResults);
   auto endCreation = std::chrono::steady_clock::now();
   auto diffCreation = endCreation - startCreation;
   std::cout << std::chrono::duration <double, std::milli> (diffCreation).count() << " ms" << std::endl;

   //Update 
   std::cout << "Updating points ... ";
   auto startUpdate = std::chrono::steady_clock::now();
   msc->updatePointsPositions(cvPoints); 
   auto endUpdate = std::chrono::steady_clock::now();
   auto diffUpdate = endUpdate - startUpdate;
   std::cout << std::chrono::duration <double, std::milli> (diffUpdate).count() << " ms" << std::endl;

   //Collide
   std::vector<bool> collideVector;
   for(size_t k=0; k<2;k++){ // 2 tet
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

   delete msc;

   std::cout << "Test : Ending.. " << std::endl << std::endl;

   return 0;
}
