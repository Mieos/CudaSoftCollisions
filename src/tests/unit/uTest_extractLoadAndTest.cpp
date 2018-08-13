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

#include <algorithm>


int main(int argc, char *argv[]){

   //Get and print the path to data
   std::string path = MIEOS_HELPERS_DATAPATH_GENERATED;

   //SavePath
   std::string savePath = MIEOS_HELPERS_DATAPATH_GENERATED;
   savePath = savePath + "/meshes/bug.ply";

   //Icosaheron
   //path = path + "/meshes/torus_vol_highres.vtk";
   path = path + "/meshes/sphere_volume.vtk";
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
   //Check results
   if(collideVector.at(1066)){
      std::cout << "BUG" << std::endl;
   } else { 
      std::cout << "OK" << std::endl;
   }

   for(size_t k=0; k<collideVector.size();k++){
      if(collideVector.at(k)){
         //std::cout << "Collide : " << k << std::endl;
      } else {
         //std::cout << "NOPE : " << k << std::endl;
      }
   }

   //DEBUG
   cv::Mat debugPoints = cv::Mat::zeros(8,3,CV_32FC1);


   size_t idTet1 = tetSelected.at(1066);
   size_t idTet2 = tetSelected.at(1109);

   vtkIdType debug_npts;
   vtkIdType* debug_pts1;
   vtkIdType* debug_pts2;
   mesh3d->GetCellPoints(idTet1, debug_npts, debug_pts1);

   std::cout << "DEBUG tet 1 = " << debug_pts1[0] << "," << debug_pts1[1] << "," << debug_pts1[2] << "," << debug_pts1[3] << std::endl;
   mesh3d->GetCellPoints(idTet2, debug_npts, debug_pts2);
   std::cout << "DEBUG tet 2 = " << debug_pts2[0] << "," << debug_pts2[1] << "," << debug_pts2[2] << "," << debug_pts2[3] << std::endl;
   
   std::cout << "TET1 = ";
   float px, py, pz;
   for(size_t k =0; k<4; k++){
      size_t idP = tetIdVector.at(1066).at(k);
      /*
      auto findT = std::find(associationResults.begin(), associationResults.end(), idTet);
      size_t idFind = findT - associationResults.begin();
      float px = pointsMat.at<float>(idFind,0);
      float py = pointsMat.at<float>(idFind,1);
      float pz = pointsMat.at<float>(idFind,2);
      */
   
      px = pointsMat.at<float>(associationResults.at(idP),0);
      py = pointsMat.at<float>(associationResults.at(idP),1);
      pz = pointsMat.at<float>(associationResults.at(idP),2);
      
      debugPoints.at<float>(k,0) = px;
      debugPoints.at<float>(k,1) = py;
      debugPoints.at<float>(k,2) = pz;
      std:: cout << associationResults.at(idP) << " , "; //"(" << idFind << ")"<< " , ";

   }
   //V2(debug)
   /*
   px = pointsMat.at<float>(associationResults.at( debug_pts1[3]),0);
   py = pointsMat.at<float>(associationResults.at( debug_pts1[3]),1);
   pz = pointsMat.at<float>(associationResults.at( debug_pts1[3]),2);
   debugPoints.at<float>(3,0) = px;
   debugPoints.at<float>(3,1) = py;
   debugPoints.at<float>(3,2) = pz;
   */ 
   std::cout << std::endl;
  

   std::cout << "TET2 = ";
   for(size_t k =0; k<4; k++){
      size_t idP = tetIdVector.at(1109).at(k);
      /*
      auto findT = std::find(associationResults.begin(), associationResults.end(), idTet);
      size_t idFind = findT - associationResults.begin();
      float px = pointsMat.at<float>(idFind,0);
      float py = pointsMat.at<float>(idFind,1);
      float pz = pointsMat.at<float>(idFind,2);
      */
      px = pointsMat.at<float>(associationResults.at(idP),0);
      py = pointsMat.at<float>(associationResults.at(idP),1);
      pz = pointsMat.at<float>(associationResults.at(idP),2);
      debugPoints.at<float>(4+k,0) = px;
      debugPoints.at<float>(4+k,1) = py;
      debugPoints.at<float>(4+k,2) = pz;

      std:: cout << associationResults.at(idP)<< "(" << idP << ")" <<  " , "; //<< "(" << idFind << ")"<< " , ";
   } 
   //V2
   /*
   px = pointsMat.at<float>(associationResults.at( debug_pts2[3]),0);
   py = pointsMat.at<float>(associationResults.at( debug_pts2[3]),1);
   pz = pointsMat.at<float>(associationResults.at( debug_pts2[3]),2);
   debugPoints.at<float>(7,0) = px;
   debugPoints.at<float>(7,1) = py;
   debugPoints.at<float>(7,2) = pz;
   */
   std::cout << std::endl;

   size_t idOrientedTetFaces[] = {
      0,1,2,
      0,2,3,
      0,3,1,
      2,1,3
   };


   float v1[3], v2[3], v3[3];

   for(size_t k=0;k<3; k++){
      v1[k] = debugPoints.at<float>(1,k) - debugPoints.at<float>(0,k);
      v2[k] = debugPoints.at<float>(2,k) - debugPoints.at<float>(0,k);
      v3[k] = debugPoints.at<float>(3,k) - debugPoints.at<float>(0,k);
   }

   double v1Vv2[3];
   v1Vv2[0] = v1[1]*v2[2] - v1[2]*v2[1];
   v1Vv2[1] = v1[2]*v2[0] - v1[0]*v2[2];
   v1Vv2[2] = v1[0]*v2[1] - v1[1]*v2[0];

   float dotP = v1Vv2[0]*v3[0] + v1Vv2[1]*v3[1] +v1Vv2[2]*v3[2];
   if(dotP>0){
      std::cout << "tet1 : WTF orientation" << std::endl;
   } else {
      std::cout << "tet1 : OK orientation" << std::endl;
   }

   for(size_t k=0;k<3; k++){
      v1[k] = debugPoints.at<float>(5,k) - debugPoints.at<float>(4,k);
      v2[k] = debugPoints.at<float>(6,k) - debugPoints.at<float>(4,k);
      v3[k] = debugPoints.at<float>(7,k) - debugPoints.at<float>(4,k);
   }

   v1Vv2[0] = v1[1]*v2[2] - v1[2]*v2[1];
   v1Vv2[1] = v1[2]*v2[0] - v1[0]*v2[2];
   v1Vv2[2] = v1[0]*v2[1] - v1[1]*v2[0];

   dotP = v1Vv2[0]*v3[0] + v1Vv2[1]*v3[1] +v1Vv2[2]*v3[2];
   if(dotP>0){
      std::cout << "tet2 : WTF orientation" << std::endl;
   } else {
      std::cout << "tet2 : OK orientation" << std::endl;
   }




   vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints > :: New();
   for(size_t k=0; k<8; k++){
      points->InsertNextPoint(debugPoints.at<float>(k,0),debugPoints.at<float>(k,1),debugPoints.at<float>(k,2));
   }
   vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
   for(size_t k=0; k<4; k++){
      vtkSmartPointer<vtkTriangle> triangleU = vtkSmartPointer<vtkTriangle>::New();
      triangleU->GetPointIds()->SetId( 0, idOrientedTetFaces[3*k] );
      triangleU->GetPointIds()->SetId( 1, idOrientedTetFaces[3*k+1] );
      triangleU->GetPointIds()->SetId( 2, idOrientedTetFaces[3*k+2] );
      triangles->InsertNextCell(triangleU);

   }
   for(size_t k=0; k<4; k++){
      vtkSmartPointer<vtkTriangle> triangleU = vtkSmartPointer<vtkTriangle>::New();
      triangleU->GetPointIds()->SetId( 0, 4+idOrientedTetFaces[3*k] );
      triangleU->GetPointIds()->SetId( 1, 4+idOrientedTetFaces[3*k+1] );
      triangleU->GetPointIds()->SetId( 2, 4+idOrientedTetFaces[3*k+2] );
      triangles->InsertNextCell(triangleU);

   }

   vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
   polyData->SetPoints(points);
   polyData->SetPolys(triangles);

   vtkSmartPointer<vtkPLYWriter> plyWriter = vtkSmartPointer<vtkPLYWriter>::New();
   plyWriter->SetFileName(savePath.c_str());
   plyWriter->SetFileTypeToASCII();

#if VTK_MAJOR_VERSION <= 5
   plyWriter->SetInput(polyData);
#else
   plyWriter->SetInputData(polyData);
#endif
   plyWriter->Write();

   delete msc;

   std::cout << "Test : Ending.. " << std::endl << std::endl;



   return 0;
}
