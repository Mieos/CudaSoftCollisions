#include <iostream>
#include <string>
#include <fstream>

#include "vtkPLYWriter.h"

#include "pathsHelpers.hpp"

#include "csc/MeshStructureExtractor.hpp"
#include "csc/MeshHelpers.hpp"

#include <map>
#include <vtkTriangle.h>
#include <vtkCellArray.h>
#include <vtkPLYWriter.h>

int main(int argc, char *argv[]){

   //Get and print the path to data
   std::string path = MIEOS_HELPERS_DATAPATH_GENERATED;

   //SavePath
   std::string savePath = MIEOS_HELPERS_DATAPATH_GENERATED;
   savePath = savePath + "/meshes/test3Dmodel.ply";

   //Icosaheron
   //path = path + "/meshes/sphere_volume.vtk";
   path = path + "/meshes/triangularDipyramid_volume.vtk";
   std::cout << "Reading: " << path << std::endl;

   //Test model extraction
   cv::Mat matPointsTet;
   std::vector<size_t> associationResults;
   std::vector<size_t> tetSelected;
   std::vector<std::vector<size_t>> tetIdVector;
   MeshStructureExtractor::extractModelFromFile(path,matPointsTet,tetIdVector,associationResults,tetSelected);

   //Read 3D mesh
   vtkSmartPointer<vtkUnstructuredGrid> msh;
   MeshHelpers::readVolumeMeshVTK(path,msh);

   bool bugDetected=false;
   //First, test if the association is correct
   for(size_t k=0; k<associationResults.size(); k++){

      float value_a[3];
      float value_b[3];
      double value_b_d[3];
      msh->GetPoints()->GetPoint(associationResults.at(k),value_b_d);

      for(size_t i=0; i<3; i++){
         value_a[i] = matPointsTet.at<float>(k,i);
         value_b[i] = float(value_b_d[i]);
         if(value_b[i] != value_a[i]){
            std::cout << "Issue: values " << value_a[i] << " and " << value_b[i] << " should be equals..." << std::endl;
            bugDetected=true;
         }
      }

   }
   if(bugDetected){
      std::cout << "Problem detected..." << std::endl;
   } else {
      std::cout << "No problem detected" << std::endl;
   }

   std::cout << "NUM points in matrix : " << matPointsTet.rows << std::endl;

   //Create the mesh
   std::map<size_t,bool> mapPoints;
   std::map<std::tuple<size_t,size_t,size_t>,std::vector<size_t>> mapFaces;
   vtkIdType npts;
   vtkIdType *pts;
   size_t idTetFaces[] = {
      0,1,2,3,
      0,1,3,2,
      0,2,3,1,
      1,2,3,0
   };
   

   
   for(size_t k=0; k<tetSelected.size(); k++){
      msh->GetCellPoints(tetSelected.at(k), npts, pts);
      for(size_t i=0; i<4; i++){
         size_t i0 = pts[idTetFaces[4*i]];
         size_t i1 = pts[idTetFaces[4*i+1]];
         size_t i2 = pts[idTetFaces[4*i+2]];
         size_t iIntruder = pts[idTetFaces[4*i+3]];
         mapPoints[i0]=true;
         mapPoints[i1]=true;
         mapPoints[i2]=true;
         std::vector<size_t> faceNums;
         faceNums.push_back(i0);
         faceNums.push_back(i1);
         faceNums.push_back(i2);
         std::sort(faceNums.begin(),faceNums.end());
         std::tuple<size_t, size_t, size_t> keyU = std::make_tuple(faceNums.at(0),faceNums.at(1),faceNums.at(2));
         if (mapFaces.find(keyU) == mapFaces.end() ) {
            std::vector<size_t> interVector;
            interVector.push_back(iIntruder);
            mapFaces[keyU]=interVector;   
         } else {
            mapFaces[keyU].push_back(iIntruder);
         }
      }
   }
   

   //Points
   vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints > :: New();
   std::vector<size_t> pointsID;
   for(std::map<size_t,bool>::iterator it = mapPoints.begin(); it != mapPoints.end(); ++it) {
      pointsID.push_back(it->first);
      double pUsedF[3];
      size_t nP = it->first;
      std::cout << "DEBUG nP = " << nP << std::endl;
      msh->GetPoints()->GetPoint(nP,pUsedF);
      if(nP==0){
         std::cout << (pUsedF[0]) << " , " << (pUsedF[1]) << " , " << (pUsedF[1]) << std::endl;
         points->InsertNextPoint(-0.0375624,0.19518,-0.132803);
      } else {
         points->InsertNextPoint(pUsedF[0],pUsedF[1],pUsedF[2]);
      }
   }

   //Faces

   vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
   for(std::map<std::tuple<size_t,size_t,size_t>,std::vector<size_t>>::iterator it = mapFaces.begin(); it != mapFaces.end(); ++it) {

      vtkSmartPointer<vtkTriangle> triangleU = vtkSmartPointer<vtkTriangle>::New();

      size_t i0 = std::get<0>(it->first);
      size_t i1 = std::get<1>(it->first);
      size_t i2 = std::get<2>(it->first);

      auto searchId0 = std::find(pointsID.begin(), pointsID.end(), i0);
      auto searchId1 = std::find(pointsID.begin(), pointsID.end(), i1);
      auto searchId2 = std::find(pointsID.begin(), pointsID.end(), i2);

      size_t i0R = searchId0 - pointsID.begin();
      size_t i1R = searchId1 - pointsID.begin();
      size_t i2R = searchId2 - pointsID.begin();

      //Check orientation
      size_t iIntruder = it->second.at(0);
      double p1Used[3], p2Used[3], p3Used[3], p4Used[3];
      msh->GetPoints()->GetPoint(i0,p1Used);
      msh->GetPoints()->GetPoint(i1,p2Used);
      msh->GetPoints()->GetPoint(i2,p3Used);
      msh->GetPoints()->GetPoint(iIntruder,p4Used);
      double v1[3],v2[3],v3[3];
      for(size_t k=0; k<3; k++){
         v1[k]=p2Used[k]-p1Used[k];
         v2[k]=p3Used[k]-p1Used[k];
         v3[k]=p4Used[k]-p1Used[k];
      }
      //Reorient
      if(MeshHelpers::checkDirectVectorOrientation(v1,v2,v3)){
         triangleU->GetPointIds()->SetId( 0, i0R );
         triangleU->GetPointIds()->SetId( 1, i2R );
         triangleU->GetPointIds()->SetId( 2, i1R );
      } else {
         triangleU->GetPointIds()->SetId( 0, i0R );
         triangleU->GetPointIds()->SetId( 1, i1R );
         triangleU->GetPointIds()->SetId( 2, i2R );
      }
      triangles->InsertNextCell(triangleU);

   }
   
   //Polydata
   vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
   polyData->SetPoints(points);
   polyData->SetPolys(triangles);

   //Save
   vtkSmartPointer<vtkPLYWriter> plyWriter = vtkSmartPointer<vtkPLYWriter>::New();
   plyWriter->SetFileName(savePath.c_str());
   //plyWriter->SetFileTypeToASCII();
   plyWriter->SetFileTypeToBinary();
#if VTK_MAJOR_VERSION <= 5
   plyWriter->SetInput(polyData);
#else
   plyWriter->SetInputData(polyData);
#endif
   plyWriter->Write();
   //

   std::cout << "Test : Ending.. " << std::endl << std::endl;



   return 0;
}
