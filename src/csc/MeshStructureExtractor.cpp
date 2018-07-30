#include "csc/MeshStructureExtractor.hpp"
#include "csc/MeshHelpers.hpp"

#include <algorithm>
#include <vtkCellArray.h>
#include <map>

//Opencv
#include "opencv2/features2d/features2d.hpp"

bool MeshStructureExtractor::extractModelFromFile(std::string fileName, cv::Mat & resultStructMatTet, std::vector<size_t> & associationVectorResult , std::vector<size_t> & tetSelected){

   vtkSmartPointer<vtkUnstructuredGrid> msh;
   if(!MeshHelpers::readVolumeMeshVTK(fileName,msh)){
      return false;
   }

   return extractModelFromMesh(msh, resultStructMatTet, associationVectorResult, tetSelected);

}

bool MeshStructureExtractor::extractModelFromMesh(const vtkSmartPointer<vtkUnstructuredGrid> & mesh3d, cv::Mat & resultStructMatTet, std::vector<size_t> & associationVectorResult, std::vector<size_t> & tetSelected){

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

   //Tetrahedron map
   cv::Mat tetrahedronMat = cv::Mat::zeros(4*associationsVect.size(),3,CV_32FC1);

   for(size_t k=0; k<associationsVect.size(); k++){
      
      size_t numTet = associationsVect.at(k).second.first;
      size_t numExtFace = associationsVect.at(k).second.second;

      tetSelected.push_back(numTet);

      //Get the corresponding tet
      mesh3d->GetCellPoints(numTet, npts, pts);
      
      //SubID of the external face
      std::vector<size_t> vectExtFace;
      vectExtFace.push_back(idTetFaces[3*numExtFace]);
      vectExtFace.push_back(idTetFaces[3*numExtFace+1]);
      vectExtFace.push_back(idTetFaces[3*numExtFace+2]);

      //Search the index of the point that is not on the surface
      size_t extPointID = 0;
      while(std::find(vectExtFace.begin(),vectExtFace.end(),extPointID)!=vectExtFace.end()){
         extPointID++;
      }
      
      //Vect
      double vect1[3];
      double vect2[3];

      //Get the real id of all points
      size_t idF1[3];
      idF1[0] = pts[idTetFaces[3*numExtFace]];
      idF1[1] = pts[idTetFaces[3*numExtFace+1]];
      idF1[2] = pts[idTetFaces[3*numExtFace+2]];
      size_t pointIDout = pts[extPointID];

      //Get the corresponding data points
      double p1Data[3];
      double p2Data[3];
      double p3Data[3];
      double pOutData[3];

      mesh3d->GetPoints()->GetPoint(idF1[0],p1Data);
      mesh3d->GetPoints()->GetPoint(idF1[1],p2Data);
      mesh3d->GetPoints()->GetPoint(idF1[2],p3Data);
      mesh3d->GetPoints()->GetPoint(pointIDout,pOutData);

      //Compute the corresponding vector (in order to check if there is no reordering to do)
      double v1[3];
      double v2[3];
      double v3[3];

      //
      for(size_t i=0;i<3; i++){
         v1[i] = p2Data[i] - p1Data[i];
         v2[i] = p3Data[i] - p1Data[i];
         v3[i] = pOutData[i] - p1Data[i];
      }

      if(MeshHelpers::checkDirectVectorOrientation(v1,v2,v3)){
     
         //OK
         for(size_t j=0; j<3; j++){
            tetrahedronMat.at<float>(4*k,j) = float(p1Data[j]);
            tetrahedronMat.at<float>(4*k+1,j) = float(p2Data[j]);
            tetrahedronMat.at<float>(4*k+2,j) = float(p3Data[j]);
            tetrahedronMat.at<float>(4*k+3,j) = float(pOutData[j]);
         }

      } else {
      
         //Reorder
         for(size_t j=0; j<3; j++){
            tetrahedronMat.at<float>(4*k,j) = float(p1Data[j]);
            tetrahedronMat.at<float>(4*k+1,j) = float(p3Data[j]);
            tetrahedronMat.at<float>(4*k+2,j) = float(p2Data[j]);
            tetrahedronMat.at<float>(4*k+3,j) = float(pOutData[j]);
         }

      }

   }

   cv::Mat matMesh3D = cv::Mat::zeros(mesh3d->GetNumberOfPoints(),3,CV_32FC1);
   
   for(size_t k=0; k<mesh3d->GetNumberOfPoints();k++){
   
      double pts3D[3];
      mesh3d->GetPoints()->GetPoint(k,pts3D);
      for(size_t i=0; i<3; i++){
         matMesh3D.at<float>(k,i) = float(pts3D[i]);
      }
   }

   //Get final associations
   cv::Ptr< cv::DescriptorMatcher > matcher = cv::DescriptorMatcher::create("BruteForce");
   std::vector< std::vector< cv::DMatch> > matches;
   matcher->knnMatch(tetrahedronMat, matMesh3D, matches, 1);
   
   for(size_t k=0; k<tetrahedronMat.rows; k++){
      size_t indexL = matches[k][0].trainIdx;
      associationVectorResult.push_back(indexL);
   }

   resultStructMatTet = tetrahedronMat.clone();

   return true;

}
