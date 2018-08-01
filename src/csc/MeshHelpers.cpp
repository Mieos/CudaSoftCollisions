#include "csc/MeshHelpers.hpp"

#include <vtkUnstructuredGridReader.h>
#include <vtkCellArray.h>

#include <map>
#include <algorithm>    // std::sort

#include <vtkTriangle.h>

//Read 3d VTK file
bool MeshHelpers::readVolumeMeshVTK(const std::string & fileName, vtkSmartPointer<vtkUnstructuredGrid> & mesh3d){

   vtkSmartPointer<vtkUnstructuredGridReader> reader3d = vtkSmartPointer<vtkUnstructuredGridReader>::New();
   reader3d->SetFileName (fileName.c_str() );
   reader3d->Update();
   mesh3d = reader3d->GetOutput();
   return true;

}

bool MeshHelpers::checkDirectVectorOrientation(const double* v1, const double* v2, const double* v3){

   //Vectoriel product
   double v1Vv2[3];
   
   v1Vv2[0] = v1[1]*v2[2] - v1[2]*v2[1];
   v1Vv2[1] = v1[2]*v2[0] - v1[0]*v2[2];
   v1Vv2[2] = v1[0]*v2[1] - v1[1]*v2[0];

   //Dot product

   double dotP = 0;
   for(size_t k=0; k<3; k++){
      dotP+= v1Vv2[k]*v3[k];
   }

   return (dotP > 0);

}

bool MeshHelpers::getSurfaceOfVolMesh(const vtkSmartPointer<vtkUnstructuredGrid> & mesh3d, vtkSmartPointer<vtkPolyData> & resultPolyData){

   //Init the surface result
   resultPolyData = vtkSmartPointer<vtkPolyData>::New();
   
   size_t numberCells = mesh3d->GetNumberOfCells(); 
   size_t numberOfPoints = mesh3d->GetNumberOfPoints();

   std::vector<float> pointsVect;

   double interPoint[3];
   for(size_t k=0; k<numberOfPoints; k++){
      mesh3d->GetPoint(k,interPoint);
      pointsVect.push_back(float(interPoint[0]));
      pointsVect.push_back(float(interPoint[1]));
      pointsVect.push_back(float(interPoint[2]));
   }

   vtkIdType npts, nfaces;
   vtkIdType *pts; 
   vtkIdType *idFacesArray;

   //std::map<std::tuple<size_t, size_t, size_t>,size_t> mapIdsFaces;
   std::map<std::tuple<size_t, size_t, size_t>, std::vector<size_t>> mapIdsFaces;
   std::map<size_t,bool> mapPoints;

   size_t idTetFaces[] = {
      0,1,2,
      0,1,3,
      0,2,3,
      1,2,3
   };

   for(size_t k=0; k<numberCells; k++){

      mesh3d->GetCellPoints(k, npts, pts);

      if(npts==4){ //only care of tetrahedrons

         //Get the 4 faces of the tetrahedron

         for(size_t i=0; i<4; i++){

            std::vector<size_t> faceIds;
            faceIds.push_back(pts[idTetFaces[3*i]]);
            faceIds.push_back(pts[idTetFaces[3*i+1]]);
            faceIds.push_back(pts[idTetFaces[3*i+2]]);
            std::sort(faceIds.begin(),faceIds.end());

            std::tuple<size_t, size_t, size_t> key_t = std::make_tuple(faceIds.at(0),faceIds.at(1),faceIds.at(2));
            if (mapIdsFaces.find(key_t) == mapIdsFaces.end() ) {
               //Update Maps
               std::vector<size_t> interVector;
               interVector.push_back(k);
               mapIdsFaces[key_t]=interVector;

            } else {
               //Update Maps
               mapIdsFaces[key_t].push_back(k);
            }

         }

      }

   }

   for(std::map<std::tuple<size_t, size_t, size_t>,std::vector<size_t>>::iterator it = mapIdsFaces.begin(); it != mapIdsFaces.end(); ++it) {
      size_t p1,p2,p3;
      p1 = std::get<0>(it->first);
      p2 = std::get<1>(it->first);
      p3 = std::get<2>(it->first);
      mapPoints[p1]=true;
      mapPoints[p2]=true;
      mapPoints[p3]=true;

   }
   
   std::vector<size_t> linkPoints; //get the corresponding ID for points of the surface 
   for(std::map<size_t,bool>::iterator it = mapPoints.begin(); it != mapPoints.end(); ++it) {
      if(it->second){
         linkPoints.push_back(it->first);
      }
   }

   //Create the mesh of the surface
   vtkSmartPointer<vtkPoints> pointsSurf =vtkSmartPointer<vtkPoints>::New();
   vtkSmartPointer<vtkCellArray> trianglesSurf = vtkSmartPointer<vtkCellArray>::New();

   for(size_t k=0; k<linkPoints.size();k++){
       
      double x = pointsVect.at(3*linkPoints.at(k));
      double y = pointsVect.at(3*linkPoints.at(k)+1);
      double z = pointsVect.at(3*linkPoints.at(k)+2);
      pointsSurf->InsertNextPoint(x, y, z);

   }
   

   
   for(std::map<std::tuple<size_t, size_t, size_t>,std::vector<size_t>>::iterator it = mapIdsFaces.begin(); it != mapIdsFaces.end(); ++it) {
      size_t p1,p2,p3;
      p1 = std::get<0>(it->first);
      p2 = std::get<1>(it->first);
      p3 = std::get<2>(it->first);
      size_t numIt = it->second.size();
      if(numIt==1){
         
         size_t idP1, idP2, idP3;

         auto resF1 = std::find(linkPoints.begin(), linkPoints.end(), p1);
         if(resF1!=linkPoints.end()){
            idP1 = resF1 - linkPoints.begin();
         } else {
            std::cout << "Malformed VTK file" << std::endl ;
            return false;
         }
        
         auto resF2 = std::find(linkPoints.begin(), linkPoints.end(), p2);
         if(resF2!=linkPoints.end()){
            idP2 = resF2 - linkPoints.begin();
         } else {
            std::cout << "Malformed VTK file" << std::endl ;
            return false;
         }

         auto resF3 = std::find(linkPoints.begin(), linkPoints.end(), p3);
         if(resF3!=linkPoints.end()){
            idP3 = resF3 - linkPoints.begin();
         } else { 
            std::cout << "Malformed VTK file" << std::endl ;
            return false;
         }

         //Check face orientation
         size_t numTetConsidered = it->second.at(0);
         vtkIdType nptsUsed;
         vtkIdType *ptsUsed;
         mesh3d->GetCellPoints(numTetConsidered, nptsUsed, ptsUsed);
         if(nptsUsed!=4){
            std::cout << "BUG surface extraction"<< std::endl;
            return false;
         }
         std::vector<size_t> vectorFaceU;
         vectorFaceU.push_back(p1);
         vectorFaceU.push_back(p2);
         vectorFaceU.push_back(p3);
         std::sort(vectorFaceU.begin(),vectorFaceU.end());
         bool foundFace=false;
         size_t numIntruder = 0;
         while(!foundFace){
           if(std::find(vectorFaceU.begin(),vectorFaceU.end(),ptsUsed[numIntruder])==vectorFaceU.end()){
              foundFace=true;
           } else {
              numIntruder++;
           } 
         }
         double p1Data[3], p2Data[3], p3Data[3], p4Data[3];
         for(size_t k=0; k<3; k++){
            p1Data[k]=double(pointsVect.at(3*p1+k));
            p2Data[k]=double(pointsVect.at(3*p2+k));
            p3Data[k]=double(pointsVect.at(3*p3+k));
            p4Data[k]=double(pointsVect.at(3*ptsUsed[numIntruder]+k));
         }
         double v1[3],v2[3],v3[3];
         for(size_t k=0; k<3; k++){
            v1[k]=p2Data[k]-p1Data[k];
            v2[k]=p3Data[k]-p1Data[k];
            v3[k]=p4Data[k]-p1Data[k];
         }

         //Add new face
         vtkSmartPointer<vtkTriangle> triangleU = vtkSmartPointer<vtkTriangle>::New();
         if(MeshHelpers::checkDirectVectorOrientation(v1,v2,v3)){ 
            triangleU->GetPointIds()->SetId ( 0, idP1 );
            triangleU->GetPointIds()->SetId ( 1, idP3 );
            triangleU->GetPointIds()->SetId ( 2, idP2 );
         } else {
            triangleU->GetPointIds()->SetId ( 0, idP1 );
            triangleU->GetPointIds()->SetId ( 1, idP2 );
            triangleU->GetPointIds()->SetId ( 2, idP3 );
         }
         trianglesSurf->InsertNextCell(triangleU);

      }
   }

   resultPolyData->SetPoints(pointsSurf);
   resultPolyData->SetPolys(trianglesSurf);

   return true;

}

