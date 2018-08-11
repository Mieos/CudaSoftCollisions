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
   std::string pathNotCollideTet1 = path + "/meshes/notIntersectingTet/tet1.ply";
   std::string pathNotCollideTet2 = path + "/meshes/notIntersectingTet/tet2.ply";

   //Intersecting
   std::string pathCollideTet1 = path + "/meshes/intersectingTet/tet1.ply"; 
   std::string pathCollideTet2 = path + "/meshes/intersectingTet/tet2.ply";

   //Not intersecting Tet1
   vtkSmartPointer<vtkPLYReader> readerNIt1 = vtkSmartPointer<vtkPLYReader>::New(); 
   readerNIt1->SetFileName ( pathNotCollideTet1.c_str() );
   readerNIt1->Update();
   vtkSmartPointer<vtkPolyData> meshNItet1 = readerNIt1->GetOutput();

   //Not intersecting Tet2 
   vtkSmartPointer<vtkPLYReader> readerNIt2 = vtkSmartPointer<vtkPLYReader>::New(); 
   readerNIt2->SetFileName ( pathNotCollideTet2.c_str() );
   readerNIt2->Update();
   vtkSmartPointer<vtkPolyData> meshNItet2 = readerNIt2->GetOutput();

   //Intersecting Tet1 
   vtkSmartPointer<vtkPLYReader> readerIt1 = vtkSmartPointer<vtkPLYReader>::New(); 
   readerIt1->SetFileName ( pathCollideTet1.c_str() );
   readerIt1->Update();
   vtkSmartPointer<vtkPolyData> meshItet1 = readerIt1->GetOutput();

   //Intersecting Tet2 
   vtkSmartPointer<vtkPLYReader> readerIt2 = vtkSmartPointer<vtkPLYReader>::New(); 
   readerIt2->SetFileName ( pathCollideTet2.c_str() );
   readerIt2->Update();
   vtkSmartPointer<vtkPolyData> meshItet2 = readerIt2->GetOutput();

   //Copy points
   size_t numPointsNItet1 = meshNItet1->GetNumberOfPoints();
   size_t numPointsNItet2 = meshNItet2->GetNumberOfPoints();
   size_t numPointsNI = numPointsNItet1 + numPointsNItet2;
   size_t numPointsItet1 = meshItet1->GetNumberOfPoints();
   size_t numPointsItet2 = meshItet2->GetNumberOfPoints();
   size_t numPointsI = numPointsItet1 + numPointsItet2;
  
   cv::Mat niPoints = cv::Mat(numPointsNI,3,CV_32FC1);
   cv::Mat iPoints = cv::Mat(numPointsI,3,CV_32FC1);
   double interPoint[3];
   //Not intersecting
   for(size_t k=0; k<meshNItet1->GetNumberOfPoints(); k++){
      meshNItet1->GetPoint(k,interPoint);
      niPoints.at<float>(k,0)=float(interPoint[0]);
      niPoints.at<float>(k,1)=float(interPoint[1]);
      niPoints.at<float>(k,2)=float(interPoint[2]);
   }
   for(size_t k=0; k<meshNItet2->GetNumberOfPoints(); k++){
      meshNItet2->GetPoint(k,interPoint);
      niPoints.at<float>(k+numPointsNItet1,0)=float(interPoint[0]);
      niPoints.at<float>(k+numPointsNItet1,1)=float(interPoint[1]);
      niPoints.at<float>(k+numPointsNItet1,2)=float(interPoint[2]);
   }
   //Intersecting
   for(size_t k=0; k<meshItet1->GetNumberOfPoints(); k++){
      meshItet1->GetPoints()->GetPoint(k,interPoint);
      iPoints.at<float>(k,0)=float(interPoint[0]);
      iPoints.at<float>(k,1)=float(interPoint[1]);
      iPoints.at<float>(k,2)=float(interPoint[2]);
   }
   for(size_t k=0; k<meshItet2->GetNumberOfPoints(); k++){
      meshItet2->GetPoints()->GetPoint(k,interPoint);
      iPoints.at<float>(k+numPointsItet1,0)=float(interPoint[0]);
      iPoints.at<float>(k+numPointsItet1,1)=float(interPoint[1]);
      iPoints.at<float>(k+numPointsItet1,2)=float(interPoint[2]);
   }

   //Create the two tet model
   float v1[3], v2[3], v1Vv2[3], v3[3];
   float dotT;
   std::vector<std::vector<size_t>> tetIdVectorI; 
   std::vector<std::vector<size_t>> tetIdVectorNI;
   //Model not intersecting
   //TET1
   for(size_t k=0; k<3; k++){
      v1[k] = niPoints.at<float>(1,k) - niPoints.at<float>(0,k);
      v2[k] = niPoints.at<float>(2,k) - niPoints.at<float>(0,k);
      v3[k] = niPoints.at<float>(3,k) - niPoints.at<float>(0,k);
   }
   v1Vv2[0] = v1[1]*v2[2] - v1[2]*v2[1];
   v1Vv2[1] = v1[2]*v2[0] - v1[0]*v2[2];
   v1Vv2[2] = v1[0]*v2[1] - v1[1]*v2[0];
   dotT = v1Vv2[0]*v3[0] + v1Vv2[1]*v3[1] + v1Vv2[2]*v3[2];
   if(dotT<0){
      std::vector<size_t> interVect;
      interVect.push_back(0);
      interVect.push_back(1);
      interVect.push_back(2);
      interVect.push_back(3);
      tetIdVectorNI.push_back(interVect);
   } else {
      std::vector<size_t> interVect;
      interVect.push_back(0);
      interVect.push_back(2);
      interVect.push_back(1);
      interVect.push_back(3);
      tetIdVectorNI.push_back(interVect);
   }
   //TET2
   for(size_t k=0; k<3; k++){
      v1[k] = niPoints.at<float>(1+numPointsNItet1,k) - niPoints.at<float>(numPointsNItet1,k);
      v2[k] = niPoints.at<float>(2+numPointsNItet1,k) - niPoints.at<float>(numPointsNItet1,k);
      v3[k] = niPoints.at<float>(3+numPointsNItet1,k) - niPoints.at<float>(numPointsNItet1,k);
   }
   v1Vv2[0] = v1[1]*v2[2] - v1[2]*v2[1];
   v1Vv2[1] = v1[2]*v2[0] - v1[0]*v2[2];
   v1Vv2[2] = v1[0]*v2[1] - v1[1]*v2[0];
   dotT = v1Vv2[0]*v3[0] + v1Vv2[1]*v3[1] + v1Vv2[2]*v3[2];
   if(dotT<0){
      std::vector<size_t> interVect;
      interVect.push_back(numPointsNItet1);
      interVect.push_back(numPointsItet1+1);
      interVect.push_back(numPointsNItet1+2);
      interVect.push_back(numPointsNItet1+3);
      tetIdVectorNI.push_back(interVect);
   } else {
      std::vector<size_t> interVect;
      interVect.push_back(numPointsNItet1);
      interVect.push_back(numPointsNItet1+2);
      interVect.push_back(numPointsNItet1+1);
      interVect.push_back(numPointsNItet1+3);
      tetIdVectorNI.push_back(interVect);
   }

   //Model intersecting
   //TET1
   for(size_t k=0; k<3; k++){
      v1[k] = iPoints.at<float>(1,k) - iPoints.at<float>(0,k);
      v2[k] = iPoints.at<float>(2,k) - iPoints.at<float>(0,k);
      v3[k] = iPoints.at<float>(3,k) - iPoints.at<float>(0,k);
   }
   v1Vv2[0] = v1[1]*v2[2] - v1[2]*v2[1];
   v1Vv2[1] = v1[2]*v2[0] - v1[0]*v2[2];
   v1Vv2[2] = v1[0]*v2[1] - v1[1]*v2[0];
   dotT = v1Vv2[0]*v3[0] + v1Vv2[1]*v3[1] + v1Vv2[2]*v3[2];
   if(dotT<0){
      std::vector<size_t> interVect;
      interVect.push_back(0);
      interVect.push_back(1);
      interVect.push_back(2);
      interVect.push_back(3);
      tetIdVectorI.push_back(interVect);
   } else {
      std::vector<size_t> interVect;
      interVect.push_back(0);
      interVect.push_back(2);
      interVect.push_back(1);
      interVect.push_back(3);
      tetIdVectorI.push_back(interVect);
   }
   //TET2
   for(size_t k=0; k<3; k++){
      v1[k] = iPoints.at<float>(1+numPointsItet1,k) - iPoints.at<float>(numPointsItet1,k);
      v2[k] = iPoints.at<float>(2+numPointsItet1,k) - iPoints.at<float>(numPointsItet1,k);
      v3[k] = iPoints.at<float>(3+numPointsItet1,k) - iPoints.at<float>(numPointsItet1,k);
   }
   v1Vv2[0] = v1[1]*v2[2] - v1[2]*v2[1];
   v1Vv2[1] = v1[2]*v2[0] - v1[0]*v2[2];
   v1Vv2[2] = v1[0]*v2[1] - v1[1]*v2[0];
   dotT = v1Vv2[0]*v3[0] + v1Vv2[1]*v3[1] + v1Vv2[2]*v3[2];
   if(dotT<0){
      std::vector<size_t> interVect;
      interVect.push_back(numPointsItet1);
      interVect.push_back(numPointsItet1+1);
      interVect.push_back(numPointsItet1+2);
      interVect.push_back(numPointsItet1+3);
      tetIdVectorI.push_back(interVect);
   } else {
      std::vector<size_t> interVect;
      interVect.push_back(numPointsItet1);
      interVect.push_back(numPointsItet1+2);
      interVect.push_back(numPointsItet1+1);
      interVect.push_back(numPointsItet1+3);
      tetIdVectorI.push_back(interVect);
   }

   //Vector association (identity)
   std::vector<size_t> associationResultsNI;
   for(size_t k=0; k<numPointsNI ; k++){
      associationResultsNI.push_back(k);
   }
   std::vector<size_t> associationResultsI;
   for(size_t k=0; k<numPointsI ; k++){ 
      associationResultsI.push_back(k);
   }

   //Create the object (NI)
   std::cout << "Creating model ... ";
   auto startCreation = std::chrono::steady_clock::now();  
   MeshStructureCollider* mscNI = new MeshStructureCollider(niPoints,tetIdVectorNI,associationResultsNI);
   auto endCreation = std::chrono::steady_clock::now();
   auto diffCreation = endCreation - startCreation;
   std::cout << std::chrono::duration <double, std::milli> (diffCreation).count() << " ms" << std::endl;

   //Update 
   std::cout << "Updating points ... ";
   auto startUpdate = std::chrono::steady_clock::now();
   mscNI->updatePointsPositions(niPoints); 
   auto endUpdate = std::chrono::steady_clock::now();
   auto diffUpdate = endUpdate - startUpdate;
   std::cout << std::chrono::duration <double, std::milli> (diffUpdate).count() << " ms" << std::endl;

   //Collide
   std::vector<bool> collideVectorNI;
   for(size_t k=0; k<2;k++){ // 2 tet
      collideVectorNI.push_back(false);
   }
   std::cout << "Colliding ... ";
   auto startCollide = std::chrono::steady_clock::now(); 
   mscNI->collide(collideVectorNI);
   auto endCollide = std::chrono::steady_clock::now();
   auto diffCollide = endCollide - startCollide;
   std::cout << std::chrono::duration <double, std::milli> (diffCollide).count() << " ms" << std::endl;

   //Check results (no intersection here)
   bool bugDetectedNI= false;
   for(size_t k=0; k<collideVectorNI.size();k++){
      if(collideVectorNI.at(k)){
         bugDetectedNI=true;
      }
   }

   //Print if there is a problem
   if(bugDetectedNI){
      std::cout << "BUG detected : intersection detected while it should not" << std::endl;
   } else {
      std::cout << "No intersection detected : OK" << std::endl;
   }

   //Create the object (I)
   std::cout << "Creating model ... ";
   auto startCreation2 = std::chrono::steady_clock::now();  
   MeshStructureCollider* mscI = new MeshStructureCollider(iPoints,tetIdVectorI,associationResultsI);
   auto endCreation2 = std::chrono::steady_clock::now();
   auto diffCreation2 = endCreation2 - startCreation2;
   std::cout << std::chrono::duration <double, std::milli> (diffCreation2).count() << " ms" << std::endl;

   //Update 
   std::cout << "Updating points ... ";
   auto startUpdate2 = std::chrono::steady_clock::now();
   mscI->updatePointsPositions(iPoints); 
   auto endUpdate2 = std::chrono::steady_clock::now();
   auto diffUpdate2 = endUpdate2 - startUpdate2;
   std::cout << std::chrono::duration <double, std::milli> (diffUpdate2).count() << " ms" << std::endl;

   //Collide
   std::vector<bool> collideVectorI;
   for(size_t k=0; k<2;k++){ // 2 tet
      collideVectorI.push_back(false);
   }
   std::cout << "Colliding ... ";
   auto startCollide2 = std::chrono::steady_clock::now(); 
   mscI->collide(collideVectorI);
   auto endCollide2 = std::chrono::steady_clock::now();
   auto diffCollide2 = endCollide2 - startCollide2;
   std::cout << std::chrono::duration <double, std::milli> (diffCollide2).count() << " ms" << std::endl;

   //Check results (intersection)
   bool bugDetectedI= false;
   for(size_t k=0; k<collideVectorI.size();k++){
      if(!collideVectorI.at(k)){
         bugDetectedI=true;
      }
   }

   //Print if there is a problem
   if(bugDetectedI){
      std::cout << "BUG detected : intersection not detected while it should have been" << std::endl;
   } else {
      std::cout << "Intersections detected : OK" << std::endl;
   }

   delete mscNI;
   delete mscI;

   std::cout << "Test : Ending.. " << std::endl << std::endl;

   return 0;
}
