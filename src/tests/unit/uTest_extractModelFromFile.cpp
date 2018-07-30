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
   path = path + "/meshes/torus_vol.vtk";
   std::cout << "Reading: " << path << std::endl;

   //Test model extraction
   cv::Mat matPointsTet;
   std::vector<size_t> associationResults;
   std::vector<size_t> tetSelected;
   MeshStructureExtractor::extractModelFromFile(path,matPointsTet,associationResults,tetSelected);

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

   //Create the mesh
   std::map<size_t,bool> mapPoints;
   std::map<std::tuple<size_t,size_t,size_t>,bool> mapFaces;
   vtkIdType npts;
   vtkIdType *pts;
   size_t idTetFaces[] = {
      0,1,2,
      0,1,3,
      0,2,3,
      1,2,3
   };

   
   for(size_t k=0; k<tetSelected.size(); k++){
      msh->GetCellPoints(tetSelected.at(k), npts, pts);
      for(size_t i=0; i<4; i++){
         size_t i0 = pts[idTetFaces[3*i]];
         size_t i1 = pts[idTetFaces[3*i+1]];
         size_t i2 = pts[idTetFaces[3*i+2]];
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
            mapFaces[keyU]=true;   
         } else {
            mapFaces[keyU]=false;
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
      msh->GetPoints()->GetPoint(nP,pUsedF);
      points->InsertNextPoint(pUsedF[0],pUsedF[1],pUsedF[2]);
   }

   //Faces

   vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
   for(std::map<std::tuple<size_t,size_t,size_t>,bool>::iterator it = mapFaces.begin(); it != mapFaces.end(); ++it) {

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

      triangleU->GetPointIds()->SetId( 0, i0R );
      triangleU->GetPointIds()->SetId( 1, i1R );
      triangleU->GetPointIds()->SetId( 2, i2R );

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
