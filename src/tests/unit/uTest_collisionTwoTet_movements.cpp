#include <iostream>
#include <string>
#include <fstream>

#include "pathsHelpers.hpp"

#include "csc/MeshStructureExtractor.hpp"
#include "csc/MeshStructureCollider.hpp"
#include "csc/MeshHelpers.hpp"

#include <vtkPLYReader.h>
#include <chrono>

#include <vtkTriangle.h>
#include <vtkPLYWriter.h>

int main(int argc, char *argv[]){

   //Get and print the path to data
   std::string path = MIEOS_HELPERS_DATAPATH_GENERATED;
   std::string savePath = MIEOS_HELPERS_DATAPATH_GENERATED;
   savePath = savePath + "/meshes/intersectingTet/debug.ply";

   //Intersecting
   std::string pathCollideTet1 = path + "/meshes/intersectingTet/tet1.ply"; 
   std::string pathCollideTet2 = path + "/meshes/intersectingTet/tet2.ply";

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
   size_t numPointsItet1 = meshItet1->GetNumberOfPoints();
   size_t numPointsItet2 = meshItet2->GetNumberOfPoints();
   size_t numPointsI = numPointsItet1 + numPointsItet2;
  
   cv::Mat iPoints = cv::Mat(numPointsI,3,CV_32FC1);
   double interPoint[3];
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
   std::vector<size_t> associationResultsI;
   for(size_t k=0; k<numPointsI ; k++){ 
      associationResultsI.push_back(k);
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
   std::vector<std::vector<float> > rectificationMovements;
   for(size_t k=0; k<2;k++){ // 2 tet
      collideVectorI.push_back(false);
      std::vector<float> interVect2;
      rectificationMovements.push_back(interVect2);
   }
   std::cout << "Colliding ... ";
   auto startCollide2 = std::chrono::steady_clock::now(); 
   mscI->collideAndGetMovements(collideVectorI,rectificationMovements);
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

   if(!bugDetectedI){
     
      vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints > :: New();
      for(size_t k=0; k<iPoints.rows; k++){
     
         size_t indexToUsed = 0;
         if(k>=4){
            indexToUsed=1;
         }
         float dx = rectificationMovements.at(indexToUsed).at(0);
         float dy = rectificationMovements.at(indexToUsed).at(1);
         float dz = rectificationMovements.at(indexToUsed).at(2);

         double x_new = double(iPoints.at<float>(k,0) + dx);
         double y_new = double(iPoints.at<float>(k,1) + dy);
         double z_new = double(iPoints.at<float>(k,2) + dz);

         points->InsertNextPoint(x_new,y_new,z_new);

      }

      size_t idTetFaces[] = {
               0,1,2,
               0,1,3,
               0,2,3,
               1,2,3
         };
      vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
      for(size_t k=0; k<2; k++){
         for(size_t i=0; i<4; i++){
            size_t i0 = tetIdVectorI.at(k).at(idTetFaces[3*i]);
            size_t i1 = tetIdVectorI.at(k).at(idTetFaces[3*i+1]);
            size_t i2 = tetIdVectorI.at(k).at(idTetFaces[3*i+2]);
            vtkSmartPointer<vtkTriangle> triangleU = vtkSmartPointer<vtkTriangle>::New();
            triangleU->GetPointIds()->SetId( 0, i0 );
            triangleU->GetPointIds()->SetId( 1, i1 );
            triangleU->GetPointIds()->SetId( 2, i2 );
            triangles->InsertNextCell(triangleU);
         }
      }

      vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
      polydata->SetPoints(points);
      polydata->SetPolys(triangles);

      vtkSmartPointer<vtkPLYWriter> plyWriter = vtkSmartPointer<vtkPLYWriter>::New();
      plyWriter->SetFileName(savePath.c_str());
      plyWriter->SetFileTypeToASCII();

#if VTK_MAJOR_VERSION <= 5
      plyWriter->SetInput(polydata);
#else
      plyWriter->SetInputData(polydata);
#endif
      plyWriter->Write();

   }

   delete mscI;

   std::cout << "Test : Ending.. " << std::endl << std::endl;

   return 0;
}
