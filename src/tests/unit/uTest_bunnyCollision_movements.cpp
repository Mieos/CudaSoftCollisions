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
#include <vtkGenericDataObjectWriter.h>

#include <vtkPointData.h>

#include <chrono>

#include <algorithm>


int main(int argc, char *argv[]){

   //Get and print the path to data
   std::string path = MIEOS_HELPERS_DATAPATH_GENERATED;

   //SavePath
   std::string savePath = MIEOS_HELPERS_DATAPATH_GENERATED;
   savePath = savePath + "/meshes/fixedIntersectingBunny.vtk";

   //Icosaheron
   path = path + "/meshes/bunny_vol.vtk";
   std::cout << "Reading: " << path << std::endl;

   //Test model extraction
   cv::Mat matPointsTet;
   std::vector<size_t> associationResults;
   std::vector<size_t> tetSelected;
   std::vector<std::vector<size_t>> tetIdVector;
   std::vector<bool> malformedTet;
   std::cout << "Extracting model ... ";
   auto startExtraction = std::chrono::steady_clock::now(); 
   MeshStructureExtractor::extractModelFromFile(path,matPointsTet,tetIdVector,associationResults,tetSelected, malformedTet);
   auto endExtraction = std::chrono::steady_clock::now();
   auto diffExtraction = endExtraction - startExtraction;
   std::cout << std::chrono::duration <double, std::milli> (diffExtraction).count() << " ms" << std::endl;

   vtkSmartPointer<vtkUnstructuredGrid> mesh3d;
   MeshHelpers::readVolumeMeshVTK(path,mesh3d);
   size_t numberOfPoints = mesh3d->GetNumberOfPoints();
   cv::Mat pointsMat = cv::Mat::zeros(numberOfPoints,3,CV_32FC1);
   double interPoint[3];
   float centerS[3];
   float upVector[0];
   centerS[0]=0.0;
   centerS[1]=0.0;
   centerS[2]=0.0;
   upVector[0]=1.0;
   upVector[1]=0.0;
   upVector[2]=0.0;
   for(size_t k=0; k<numberOfPoints; k++){
      mesh3d->GetPoint(k,interPoint);
      pointsMat.at<float>(k,0)=float(interPoint[0]);
      pointsMat.at<float>(k,1)=float(interPoint[1]);
      pointsMat.at<float>(k,2)=float(interPoint[2]);
      centerS[0]+=pointsMat.at<float>(k,0);
      centerS[1]+=pointsMat.at<float>(k,0);
      centerS[2]+=pointsMat.at<float>(k,0);
   }
   centerS[0]=centerS[0]/numberOfPoints;
   centerS[1]=centerS[1]/numberOfPoints;
   centerS[2]=centerS[2]/numberOfPoints;

   float collisionNorm=0.03;
   cv::Mat pointsMatUpdateCollision = cv::Mat::zeros(numberOfPoints,3,CV_32FC1);
   float vT[3];
   float dotTest;
   for(size_t k=0; k<numberOfPoints; k++){
      vT[0] = pointsMat.at<float>(k,0) - centerS[0];
      vT[1] = pointsMat.at<float>(k,1) - centerS[0];
      vT[2] = pointsMat.at<float>(k,2) - centerS[0];
      dotTest = vT[0]*upVector[0] + vT[1]*upVector[1] + vT[2]*upVector[2];
      if(dotTest<0){
         pointsMatUpdateCollision.at<float>(k,0)=pointsMat.at<float>(k,0);
         pointsMatUpdateCollision.at<float>(k,1)=pointsMat.at<float>(k,1);
         pointsMatUpdateCollision.at<float>(k,2)=pointsMat.at<float>(k,2);
      } else {
         pointsMatUpdateCollision.at<float>(k,0)=pointsMat.at<float>(k,0)-collisionNorm*upVector[0];
         pointsMatUpdateCollision.at<float>(k,1)=pointsMat.at<float>(k,1)-collisionNorm*upVector[1];
         pointsMatUpdateCollision.at<float>(k,2)=pointsMat.at<float>(k,2)-collisionNorm*upVector[2];
      }
   }
   //Create the object
   std::cout << "Creating model ... ";
   auto startCreation = std::chrono::steady_clock::now();  
   MeshStructureCollider* msc = new MeshStructureCollider(matPointsTet,tetIdVector,associationResults,malformedTet);
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
   std::vector<std::vector<float> > rectificationMovements;
   for(size_t k=0; k<tetSelected.size();k++){
      collideVector.push_back(false);
      std::vector<float> interVect2;
      rectificationMovements.push_back(interVect2);
   }

   std::cout << "Colliding ... ";
   auto startCollide = std::chrono::steady_clock::now(); 
   msc->collide(collideVector);
   auto endCollide = std::chrono::steady_clock::now();
   auto diffCollide = endCollide - startCollide;
   std::cout << std::chrono::duration <double, std::milli> (diffCollide).count() << " ms" << std::endl;

   //Check results
   bool bugDetected=false;
   for(size_t k=0; k<collideVector.size();k++){
      if(collideVector.at(k)){
         bugDetected=true;
      }
   }

   if(bugDetected){
      std::cout << "BUG detected : intersection detected while it should not" << std::endl;
   } else {
      std::cout << "No intersection detected : OK" << std::endl;
   }

   //Add collision
   bugDetected=true;

   //Update 
   std::cout << "Updating points ... ";
   auto startUpdate2 = std::chrono::steady_clock::now();
   msc->updatePointsPositions(pointsMatUpdateCollision); 
   auto endUpdate2 = std::chrono::steady_clock::now();
   auto diffUpdate2 = endUpdate2 - startUpdate2;
   std::cout << std::chrono::duration <double, std::milli> (diffUpdate2).count() << " ms" << std::endl;

   std::cout << "Colliding ... ";
   auto startCollide2 = std::chrono::steady_clock::now(); 
   msc->collideAndGetMovements(collideVector,rectificationMovements);
   auto endCollide2 = std::chrono::steady_clock::now();
   auto diffCollide2 = endCollide2 - startCollide2;
   std::cout << std::chrono::duration <double, std::milli> (diffCollide2).count() << " ms" << std::endl;

   for(size_t k=0; k<collideVector.size();k++){
      if(collideVector.at(k)){
         bugDetected=false;
      }
   }

   if(bugDetected){
      std::cout << "Something went wrong, no intersection detected" << std::endl;
   } else { 
      std::cout << "Intersections detected : OK" << std::endl;
   }

   bool fixedMesh=false;
   size_t nbIt=0;
   size_t maxIt = 10;

   while(!fixedMesh){

      nbIt++;
      std::cout << nbIt << " / " << maxIt << std::endl;

      for(size_t k=0; k<collideVector.size(); k++){
         if(collideVector.at(k)){
            float dx = rectificationMovements.at(k).at(0);
            float dy = rectificationMovements.at(k).at(1);
            float dz = rectificationMovements.at(k).at(2);
            //std::cout << "dx = " << dx <<std::endl;
            for(size_t i=0;i<4;i++){
               size_t idToChange = associationResults.at(tetIdVector.at(k).at(i));
               pointsMatUpdateCollision.at<float>(idToChange,0)=pointsMatUpdateCollision.at<float>(idToChange,0)+dx;
               pointsMatUpdateCollision.at<float>(idToChange,1)=pointsMatUpdateCollision.at<float>(idToChange,1)+dy;
               pointsMatUpdateCollision.at<float>(idToChange,2)=pointsMatUpdateCollision.at<float>(idToChange,2)+dz;
            }
         }
      }

      //Update
      msc->updatePointsPositions(pointsMatUpdateCollision);

      //Collide
      msc->collideAndGetMovements(collideVector,rectificationMovements);

      if(nbIt>=maxIt){
         fixedMesh=true;
      }

   }

   //

   vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints > :: New();
   double x_u, y_u, z_u;
   std::vector<bool> collidingPoints;
   for(size_t k=0; k<associationResults.size(); k++){
      x_u=double(pointsMatUpdateCollision.at<float>(associationResults.at(k),0));
      y_u=double(pointsMatUpdateCollision.at<float>(associationResults.at(k),1));
      z_u=double(pointsMatUpdateCollision.at<float>(associationResults.at(k),2));
      points->InsertNextPoint(x_u,y_u,z_u);
      collidingPoints.push_back(false);
   }

   size_t idTetFaces[] = {
      0,1,2,
      0,1,3,
      0,2,3,
      1,2,3
   };
   vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
   for(size_t k=0; k<tetIdVector.size(); k++){
      bool isColliding = collideVector.at(k);
      for(size_t i=0; i<4; i++){
         size_t i0 = tetIdVector.at(k).at(idTetFaces[3*i]);
         size_t i1 = tetIdVector.at(k).at(idTetFaces[3*i+1]);
         size_t i2 = tetIdVector.at(k).at(idTetFaces[3*i+2]);
         vtkSmartPointer<vtkTriangle> triangleU = vtkSmartPointer<vtkTriangle>::New();
         //New face
         triangleU->GetPointIds()->SetId( 0, i0 );
         triangleU->GetPointIds()->SetId( 1, i1 );
         triangleU->GetPointIds()->SetId( 2, i2 );
         if(isColliding){
            collidingPoints.at(i0)=true;
            collidingPoints.at(i1)=true;
            collidingPoints.at(i2)=true;
         }
         triangles->InsertNextCell(triangleU);
      }
   }

   vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
   colors->SetNumberOfComponents(3);
   colors->SetName ("Colors");

   unsigned char red[3] = {255, 0, 0};
   unsigned char green[3] = {0, 255, 0};

   for(size_t k=0; k<collidingPoints.size(); k++){
      if(collidingPoints.at(k)){
#if VTK_MAJOR_VERSION < 7      
         colors->InsertNextTupleValue(red);
#else
         colors->InsertNextTypedTuple(red);
#endif
      } else {

#if VTK_MAJOR_VERSION < 7      
         colors->InsertNextTupleValue(green);
#else
         colors->InsertNextTypedTuple(green);
#endif
      }
   }

   vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
   polydata->SetPoints(points);
   polydata->SetPolys(triangles);
   polydata->GetPointData()->SetScalars(colors);

   vtkSmartPointer<vtkGenericDataObjectWriter> vtkWriter = vtkSmartPointer<vtkGenericDataObjectWriter>::New();
   vtkWriter->SetFileName(savePath.c_str());
   vtkWriter->SetFileTypeToASCII();

#if VTK_MAJOR_VERSION <= 5
   vtkWriter->SetInput(polydata);
#else
   vtkWriter->SetInputData(polydata);
#endif
   vtkWriter->Write();

   delete msc;

   std::cout << "Test : Ending.. " << std::endl << std::endl;


   return 0;
}
