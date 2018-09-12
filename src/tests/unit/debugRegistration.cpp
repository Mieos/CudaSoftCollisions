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

#include <vtkGenericDataObjectWriter.h>

#include <vtkPointData.h>

#include <chrono>

#include <algorithm>

#include "opencv2/features2d/features2d.hpp"


int main(int argc, char *argv[]){

   //Get and print the path to data
   std::string path = MIEOS_HELPERS_DATAPATH_GENERATED;

   //SavePath
   std::string savePath = "/media/rmodrzejewski/Maxtor/data/intersections/debug/debugR.vtk";

   //Icosaheron
   path = "/media/rmodrzejewski/Maxtor/data/metalPigsXP/model/segmentations/volume.vtk";
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
   
   for(size_t k=0; k<numberOfPoints; k++){
      mesh3d->GetPoint(k,interPoint);
      pointsMat.at<float>(k,0)=float(interPoint[0]);
      pointsMat.at<float>(k,1)=float(interPoint[1]);
      pointsMat.at<float>(k,2)=float(interPoint[2]);
   }

   //Create the object
   std::cout << "Creating model ... ";
   auto startCreation = std::chrono::steady_clock::now();  
   MeshStructureCollider* msc = new MeshStructureCollider(matPointsTet,tetIdVector,associationResults, malformedTet);
   auto endCreation = std::chrono::steady_clock::now();
   auto diffCreation = endCreation - startCreation;
   std::cout << std::chrono::duration <double, std::milli> (diffCreation).count() << " ms" << std::endl;

   cv::FileStorage fs;
   cv::Mat M;
   fs.open("/home/rmodrzejewski/Projects/mRegistrationDemo/build/debug.yaml", cv::FileStorage::READ);
   fs["M"] >> M;

   //Update 
   std::cout << "Updating points ... ";
   auto startUpdate = std::chrono::steady_clock::now();
   msc->updatePointsPositions(M); 
   auto endUpdate = std::chrono::steady_clock::now();
   auto diffUpdate = endUpdate - startUpdate;
   std::cout << std::chrono::duration <double, std::milli> (diffUpdate).count() << " ms" << std::endl;

   //Collide
   std::vector<bool> collideVector;
   std::vector<bool> invertVector; 
   std::vector<std::vector<float> >  movVect1;
   std::vector<std::vector<float> >  movVect2;

   for(size_t k=0; k<tetSelected.size();k++){
      
      std::vector<float> intervect1;
      std::vector<float> intervect2;
      movVect1.push_back(intervect1);
      movVect2.push_back(intervect2);

      collideVector.push_back(false);
      invertVector.push_back(false);

   }
   
   std::cout << "Colliding ... ";
   auto startCollide = std::chrono::steady_clock::now(); 
   //msc->collide(collideVector);
   msc->collideAndGetMovementsAndInversions(collideVector, movVect1, invertVector, movVect2);
   auto endCollide = std::chrono::steady_clock::now();
   auto diffCollide = endCollide - startCollide;
   std::cout << std::chrono::duration <double, std::milli> (diffCollide).count() << " ms" << std::endl;
 
   //Tet faces
   size_t idTetFaces[] = {
         0,1,2,3,
         0,1,3,2,
         0,2,3,1,
         1,2,3,0
   };

   vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints > :: New();
   double x_u, y_u, z_u;
   
   //Result vect
   std::vector<bool> collidingPoints;
   std::vector<bool> invertingPoints;

   for(size_t k=0; k<associationResults.size(); k++){
      
      x_u=double(M.at<cv::Point3d>(associationResults.at(k)).x);
      y_u=double(M.at<cv::Point3d>(associationResults.at(k)).y);
      z_u=double(M.at<cv::Point3d>(associationResults.at(k)).z);
      
      points->InsertNextPoint(x_u,y_u,z_u);
      
      collidingPoints.push_back(false);
      invertingPoints.push_back(false);

   }

   vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
   for(size_t k=0; k<tetIdVector.size(); k++){
      bool isColliding = collideVector.at(k);
      bool isInverted = invertVector.at(k);
      for(size_t i=0; i<4; i++){
         size_t i0 = tetIdVector.at(k).at(idTetFaces[4*i]);
         size_t i1 = tetIdVector.at(k).at(idTetFaces[4*i+1]);
         size_t i2 = tetIdVector.at(k).at(idTetFaces[4*i+2]);
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
         if(isInverted){
            invertingPoints.at(i0)=true;
            invertingPoints.at(i1)=true;
            invertingPoints.at(i2)=true;
         }
         triangles->InsertNextCell(triangleU);
      }
   }

   vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
   colors->SetNumberOfComponents(3);
   colors->SetName ("Colors");

   unsigned char green[3] = {0, 255, 0};
   unsigned char red[3] = {255, 0, 0};
   unsigned char blue[3] = {0, 0, 255};
   unsigned char purple[3] = {255, 0, 255};

   for(size_t k=0; k<collidingPoints.size(); k++){
      
      //Green
      if(!collidingPoints.at(k) && !invertingPoints.at(k)){

#if VTK_MAJOR_VERSION < 7      
         colors->InsertNextTupleValue(green);
#else
         colors->InsertNextTypedTuple(green);
#endif

      } else if(collidingPoints.at(k) && !invertingPoints.at(k)){ //Red

#if VTK_MAJOR_VERSION < 7      
         colors->InsertNextTupleValue(red);
#else
         colors->InsertNextTypedTuple(red);
#endif

      } else if(!collidingPoints.at(k) && invertingPoints.at(k)){ //Blue
 
#if VTK_MAJOR_VERSION < 7      
         colors->InsertNextTupleValue(blue);
#else
         colors->InsertNextTypedTuple(blue);
#endif
     
      } else if(collidingPoints.at(k) && invertingPoints.at(k)){
       
#if VTK_MAJOR_VERSION < 7      
         colors->InsertNextTupleValue(purple);
#else
         colors->InsertNextTypedTuple(purple);
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
