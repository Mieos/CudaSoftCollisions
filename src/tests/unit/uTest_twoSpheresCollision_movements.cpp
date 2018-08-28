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
   savePath = savePath + "/meshes/twoSpheresMovements.vtk";

   //Icosaheron
   std::string pathS1 = path + "/meshes/intersectingSpheres/s1.vtk";
   std::string pathS2 = path + "/meshes/intersectingSpheres/s2.vtk";

   //Test model extraction
   cv::Mat matPointsTet1, matPointsTet2, matPointsTetT;
   std::vector<size_t> associationResults1, associationResults2, associationResultsT;
   std::vector<size_t> tetSelected1, tetSelected2;
   std::vector<std::vector<size_t>> tetIdVector1, tetIdVector2, tetIdVectorT;
   std::cout << "Extracting model 1 ... ";
   MeshStructureExtractor::extractModelFromFile(pathS1,matPointsTet1,tetIdVector1,associationResults1,tetSelected1);
   std::cout << std::endl;
   std::cout << "Extracting model 2 ... ";
   MeshStructureExtractor::extractModelFromFile(pathS2,matPointsTet2,tetIdVector2,associationResults2,tetSelected2);
   std::cout << std::endl;

   //Meshes
   vtkSmartPointer<vtkUnstructuredGrid> mesh3d1, mesh3d2;
   MeshHelpers::readVolumeMeshVTK(pathS1,mesh3d1);
   MeshHelpers::readVolumeMeshVTK(pathS2,mesh3d2);

   //cv::Mat
   size_t numberOfPoints1 = mesh3d1->GetNumberOfPoints();
   size_t numberOfPoints2 = mesh3d2->GetNumberOfPoints();
  
   std::cout << "Number of points in mesh 1 : " << numberOfPoints1 << std::endl;
   std::cout << "Number of points in mesh 2 : " << numberOfPoints2 << std::endl;
   
   cv::Mat pointsMat = cv::Mat::zeros(numberOfPoints1+numberOfPoints2,3,CV_32FC1);
   double interPoint[3];
   //Add Mesh 1
   for(size_t k=0; k<numberOfPoints1; k++){
      mesh3d1->GetPoint(k,interPoint);
      pointsMat.at<float>(k,0)=float(interPoint[0]);
      pointsMat.at<float>(k,1)=float(interPoint[1]);
      pointsMat.at<float>(k,2)=float(interPoint[2]);
   }
   //Add Mesh 2
   for(size_t k=0; k<numberOfPoints2; k++){
      mesh3d2->GetPoint(k,interPoint);
      pointsMat.at<float>(k+numberOfPoints1,0)=float(interPoint[0]);
      pointsMat.at<float>(k+numberOfPoints1,1)=float(interPoint[1]);
      pointsMat.at<float>(k+numberOfPoints1,2)=float(interPoint[2]);
   }

   //Get mat points
   size_t numberOfPointsInBuffer1 = matPointsTet1.rows;
   size_t numberOfPointsInBuffer2 = matPointsTet2.rows;
   std::cout << "Number of points in first buffer : " << numberOfPointsInBuffer1 << std::endl;
   std::cout << "Number of points in second buffer : " << numberOfPointsInBuffer2 << std::endl;
   vconcat(matPointsTet1, matPointsTet2, matPointsTetT);

   //Get tetID vector
   for(size_t k=0; k<tetIdVector1.size(); k++){
      std::vector<size_t> vectInter;
      vectInter.push_back(tetIdVector1.at(k).at(0));
      vectInter.push_back(tetIdVector1.at(k).at(1));
      vectInter.push_back(tetIdVector1.at(k).at(2));
      vectInter.push_back(tetIdVector1.at(k).at(3));
      tetIdVectorT.push_back(vectInter);
   }
   for(size_t k=0; k<tetIdVector2.size(); k++){
      std::vector<size_t> vectInter;
      vectInter.push_back(tetIdVector2.at(k).at(0)+numberOfPointsInBuffer1);
      vectInter.push_back(tetIdVector2.at(k).at(1)+numberOfPointsInBuffer1);
      vectInter.push_back(tetIdVector2.at(k).at(2)+numberOfPointsInBuffer1);
      vectInter.push_back(tetIdVector2.at(k).at(3)+numberOfPointsInBuffer1);
      tetIdVectorT.push_back(vectInter);
   }
 
   //Get association vector
   for(size_t k=0; k<associationResults1.size(); k++){
      associationResultsT.push_back(associationResults1.at(k));
   }
   for(size_t k=0; k<associationResults2.size(); k++){
      associationResultsT.push_back(associationResults2.at(k)+numberOfPoints1);
   }

   //Create the object
   std::cout << "Creating model ... ";
   MeshStructureCollider* msc = new MeshStructureCollider(matPointsTetT,tetIdVectorT,associationResultsT);
   std::cout << std::endl;

   //Update 
   std::cout << "Updating points ... ";
   msc->updatePointsPositions(pointsMat); 
   std::cout << std::endl;

   //Collide
   std::vector<bool> collideVector;
   std::vector<std::vector<float> > rectificationMovements;
   for(size_t k=0; k<tetIdVectorT.size();k++){
      collideVector.push_back(false);
      std::vector<float> interVect2;
      rectificationMovements.push_back(interVect2);
   }

   std::cout << "Colliding ... ";
   auto startCollide = std::chrono::steady_clock::now(); 
   msc->collideAndGetMovements(collideVector,rectificationMovements);
   auto endCollide = std::chrono::steady_clock::now();
   auto diffCollide = endCollide - startCollide;
   std::cout << std::chrono::duration <double, std::milli> (diffCollide).count() << " ms" << std::endl;

   //Check results
   bool bugDetected=true;
   for(size_t k=0; k<collideVector.size();k++){
      if(collideVector.at(k)){
         bugDetected=false;
      }
   }

   
   if(bugDetected){
      std::cout << "BUG detected : no intersection detected" << std::endl;
   } else {
      std::cout << "Intersection detected : OK" << std::endl;
   }

   bool fixedMesh=false;
   size_t nbIt=0;
   size_t maxIt = 10;

   while(!fixedMesh){

      cv::Mat modificationMat = cv::Mat::zeros(pointsMat.rows,pointsMat.cols, CV_32FC1);
      cv::Mat modificationMatNumber = cv::Mat::zeros(pointsMat.rows,1, CV_8UC1);
   
      nbIt++;
      std::cout << nbIt << " / " << maxIt << std::endl;

      for(size_t k=0; k<collideVector.size(); k++){
         if(collideVector.at(k)){
            float dx = rectificationMovements.at(k).at(0);
            float dy = rectificationMovements.at(k).at(1);
            float dz = rectificationMovements.at(k).at(2);
            for(size_t i=0;i<4;i++){
               size_t idToChange = associationResultsT.at(tetIdVectorT.at(k).at(i));
               modificationMat.at<float>(idToChange,0)=modificationMat.at<float>(idToChange,0)+dx;
               modificationMat.at<float>(idToChange,1)=modificationMat.at<float>(idToChange,1)+dy;
               modificationMat.at<float>(idToChange,2)=modificationMat.at<float>(idToChange,2)+dz;
               modificationMatNumber.at<uchar>(idToChange,0) = modificationMatNumber.at<uchar>(idToChange,0)+1;
            }
         }
      }

      for(size_t k=0; k<modificationMatNumber.rows; k++){
         if(modificationMatNumber.at<uchar>(k,0)>0){
            pointsMat.at<float>(k,0)=pointsMat.at<float>(k,0)+(modificationMat.at<float>(k,0)/float(modificationMatNumber.at<uchar>(k,0)));
            pointsMat.at<float>(k,1)=pointsMat.at<float>(k,1)+(modificationMat.at<float>(k,1)/float(modificationMatNumber.at<uchar>(k,0)));
            pointsMat.at<float>(k,2)=pointsMat.at<float>(k,2)+(modificationMat.at<float>(k,2)/float(modificationMatNumber.at<uchar>(k,0)));
         }
      }

      //Update
      msc->updatePointsPositions(pointsMat);

      //Collide
      msc->collideAndGetMovements(collideVector,rectificationMovements);

      bool loopStop=true;
      for(size_t k=0; k<collideVector.size(); k++){
         if(collideVector.at(k)){
            loopStop=false;
         }
      }

      if(loopStop){
         fixedMesh=true;
      }

      if(nbIt>=maxIt){
         fixedMesh=true;
      }

   }

   //

   vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints > :: New();
   double x_u, y_u, z_u;
   std::vector<bool> collidingPoints;
   for(size_t k=0; k<associationResultsT.size(); k++){
      x_u=double(pointsMat.at<float>(associationResultsT.at(k),0));
      y_u=double(pointsMat.at<float>(associationResultsT.at(k),1));
      z_u=double(pointsMat.at<float>(associationResultsT.at(k),2));
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
   for(size_t k=0; k<tetIdVectorT.size(); k++){
      bool isColliding = collideVector.at(k);
      for(size_t i=0; i<4; i++){
         size_t i0 = tetIdVectorT.at(k).at(idTetFaces[3*i]);
         size_t i1 = tetIdVectorT.at(k).at(idTetFaces[3*i+1]);
         size_t i2 = tetIdVectorT.at(k).at(idTetFaces[3*i+2]);
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
