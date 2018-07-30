#include "csc/MeshStructureExtractor.hpp"
#include "csc/MeshHelpers.hpp"

#include <algorithm>
#include <vtkCellArray.h>
#include <map>

bool MeshStructureExtractor::extractModelFromFile(std::string fileName){

   vtkSmartPointer<vtkUnstructuredGrid> msh;
   if(!MeshHelpers::readVolumeMeshVTK(fileName,msh)){
      return false;
   }

   return extractModelFromMesh(msh);

}

bool MeshStructureExtractor::extractModelFromMesh(const vtkSmartPointer<vtkUnstructuredGrid> & mesh3d){

   //First get the surface of the model
   vtkSmartPointer<vtkPolyData> surfPolyData;
   MeshHelpers::getSurfaceOfVolMesh(mesh3d, surfPolyData);

   size_t numFaces = surfPolyData->GetNumberOfPolys();
   size_t numberCells = mesh3d->GetNumberOfCells();

   std::map<std::tuple<size_t, size_t, size_t>,size_t> mapIdsFaces;
   vtkIdType npts;
   vtkIdType *pts;

   int h; 

   //Extract the faces
   for(size_t k=0; k<numFaces; k++){

      h=surfPolyData->GetPolys()->GetNextCell(npts,pts); 

      if(h==0){
         std::cout << "Malformed polydata: escaping" << std::endl;
         return false;
      }

      if(npts==3){
         std::vector<size_t> interVector;
         interVector.push_back(pts[0]);
         interVector.push_back(pts[1]);
         interVector.push_back(pts[2]);
         std::sort(interVector.begin(),interVector.end());
         std::tuple<size_t, size_t, size_t> keyU = std::make_tuple(interVector.at(0),interVector.at(1),interVector.at(2));
         mapIdsFaces[keyU]=k;
      } else {
         std::cout << "Malformed polydata: escaping" << std::endl;
         return false;
      } 

   }

   size_t idTetFaces[] = {
      0,1,2,
      0,1,3,
      0,2,3,
      1,2,3
   };

   std::vector<std::pair<size_t,std::pair<size_t,size_t>>> associationsVect;

   //Get the associations
   for(size_t k=0; k<numberCells; k++){

      mesh3d->GetCellPoints(k, npts, pts);
      
      if(npts==4){ //only care of tetrahedrons

         for(size_t i=0; i<4; i++){

            std::vector<size_t> faceIds;
            faceIds.push_back(pts[idTetFaces[3*i]]);
            faceIds.push_back(pts[idTetFaces[3*i+1]]);
            faceIds.push_back(pts[idTetFaces[3*i+2]]);
            std::sort(faceIds.begin(),faceIds.end());

            std::tuple<size_t, size_t, size_t> key_t = std::make_tuple(faceIds.at(0),faceIds.at(1),faceIds.at(2));
            
            if (mapIdsFaces.find(key_t) == mapIdsFaces.end() ) {
               std::pair<size_t,size_t> tetPart = std::make_pair(k,i);
               std::pair<size_t,std::pair<size_t,size_t>> association = std::make_pair(mapIdsFaces[key_t],tetPart);
               associationsVect.push_back(association);
            }

         }

      }

   } 

   std::map<size_t,bool> mapPointsUsed;

   for(size_t k=0; k<associationsVect.size(); k++){
      
      size_t numTet = associationsVect.at(k).second.first;
      size_t numExtFace = associationsVect.at(k).second.second;

      mesh3d->GetCellPoints(numTet, npts, pts);
      
      std::vector<size_t> vectExtFace;
      vectExtFace.push_back(idTetFaces[3*numExtFace]);
      vectExtFace.push_back(idTetFaces[3*numExtFace+1]);
      vectExtFace.push_back(idTetFaces[3*numExtFace+2]);

      size_t otherPoints=0;
      while(std::find(vectExtFace.begin(),vectExtFace.end(),otherPoints)==vectExtFace.end()){
         otherPoints++;
      }

      //Vect
      size_t vect1[3];
      size_t vect2[3];

      

      std::cout << "F1 = "  << idTetFaces[3*numExtFace] << " , " <<
         idTetFaces[3*numExtFace+1] << " , " <<
         idTetFaces[3*numExtFace+2] << std::endl;
      std::cout << "Debug = " << otherPoints << std::endl;

   }

   return false; //TODO

}
