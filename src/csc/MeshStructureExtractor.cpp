#include "csc/MeshStructureExtractor.hpp"
#include "csc/MeshHelpers.hpp"

#include <algorithm>
#include <vtkCellArray.h>
#include <map>

//Opencv
#include "opencv2/features2d/features2d.hpp"

bool MeshStructureExtractor::extractModelFromFile(std::string fileName, cv::Mat & resultStructMatTet, std::vector<std::vector<size_t>> & tetIdVector, std::vector<size_t> & associationVectorResult , std::vector<size_t> & tetSelected){

   vtkSmartPointer<vtkUnstructuredGrid> msh;
   if(!MeshHelpers::readVolumeMeshVTK(fileName,msh)){
      return false;
   }

   return extractModelFromMesh(msh, resultStructMatTet, tetIdVector, associationVectorResult, tetSelected);

}

bool MeshStructureExtractor::extractModelFromMesh(const vtkSmartPointer<vtkUnstructuredGrid> & mesh3d, cv::Mat & resultStructMatTet, std::vector<std::vector<size_t>> & tetIdVector, std::vector<size_t> & associationVectorResult, std::vector<size_t> & tetSelected){

   //First get the surface of the model
   vtkSmartPointer<vtkPolyData> surfPolyData;
   MeshHelpers::getSurfaceOfVolMesh(mesh3d, surfPolyData);

   size_t numFaces = surfPolyData->GetNumberOfPolys();
   size_t numberCells = mesh3d->GetNumberOfCells();

   vtkIdType npts;
   vtkIdType *pts;

   int h; 

   std::vector<std::tuple<size_t, size_t, size_t>> vectorIdFacesSurface;
   //Extract the faces of the surfaces
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
         vectorIdFacesSurface.push_back(keyU);
         //mapIdsFaces[keyU]=k;
      } else {
         std::cout << "Malformed polydata: escaping" << std::endl;
         return false;
      } 

   }

   //Extract the points of the surface
   cv::Mat matSurfPoints3D = cv::Mat::zeros(surfPolyData->GetNumberOfPoints(),3,CV_32FC1);
   for(size_t k = 0; k<surfPolyData->GetNumberOfPoints(); k++){
      double pInter[3];
      surfPolyData->GetPoint(k,pInter);
      for(size_t i=0; i<3; i++){
         matSurfPoints3D.at<float>(k,i) = float(pInter[i]);
      }
   }

   //Extract face of the volumetric mesh
   size_t idTetFaces[] = {
      0,1,2,
      0,1,3,
      0,2,3,
      1,2,3
   };
   std::map<std::tuple<size_t, size_t, size_t>,std::vector<size_t>> mapIdsFacesVolum;
   for(size_t k=0; k<numberCells; k++){
      mesh3d->GetCellPoints(k, npts, pts);
      if(npts==4){ //only care of tetrahedrons
         for(size_t i=0; i<4; i++){
            std::vector<size_t> faceIds;
            faceIds.push_back(pts[idTetFaces[3*i]]);
            faceIds.push_back(pts[idTetFaces[3*i+1]]);
            faceIds.push_back(pts[idTetFaces[3*i+2]]);
            std::sort(faceIds.begin(),faceIds.end());
            std::tuple<size_t, size_t, size_t> keyU = std::make_tuple(faceIds.at(0), faceIds.at(1), faceIds.at(2));
            if (mapIdsFacesVolum.find(keyU) == mapIdsFacesVolum.end() ) {
               std::vector<size_t> interVector;
               interVector.push_back(k);
               mapIdsFacesVolum[keyU]=interVector;
            } else {
               mapIdsFacesVolum[keyU].push_back(k);
            }
         }
      }
   }

   //Extract the points of the volumetric data
   cv::Mat matVolumePoints3D = cv::Mat::zeros(mesh3d->GetNumberOfPoints(),3,CV_32FC1);
   for(size_t k=0; k<mesh3d->GetNumberOfPoints(); k++){
      double pInter[3];
      mesh3d->GetPoint(k,pInter);
      for(size_t i=0; i<3; i++){
         matVolumePoints3D.at<float>(k,i) = float(pInter[i]);
      }
   }

   //Get the association surface/Volume
   std::vector<size_t> associationsSV;
   cv::Ptr< cv::DescriptorMatcher > matcherSV = cv::DescriptorMatcher::create("BruteForce");
   std::vector< std::vector< cv::DMatch> > matchesSV;
   matcherSV->knnMatch(matSurfPoints3D, matVolumePoints3D, matchesSV, 1);
   for(size_t k=0; k<matSurfPoints3D.rows; k++){
      size_t indexL = matchesSV[k][0].trainIdx;
      associationsSV.push_back(indexL);
   }

   //Search the faces of the surface
   std::vector<std::pair<size_t,size_t>> assoFaceSFaceV;
   for(size_t k=0; k<vectorIdFacesSurface.size(); k++){
      
      //Get the face
      size_t p1 = std::get<0>(vectorIdFacesSurface.at(k));
      size_t p2 = std::get<1>(vectorIdFacesSurface.at(k));
      size_t p3 = std::get<2>(vectorIdFacesSurface.at(k));
      //Get the corresponding point in volumetric reference
      auto p1C = associationsSV.at(p1);
      auto p2C = associationsSV.at(p2);
      auto p3C = associationsSV.at(p3);
      std::vector<size_t> faceTranslated;
      faceTranslated.push_back(p1C);
      faceTranslated.push_back(p2C);
      faceTranslated.push_back(p3C);
      std::sort(faceTranslated.begin(),faceTranslated.end());
      std::tuple<size_t, size_t, size_t> keyU = std::make_tuple(faceTranslated.at(0), faceTranslated.at(1), faceTranslated.at(2));

      //Search the face
      if (mapIdsFacesVolum.find(keyU) == mapIdsFacesVolum.end() ) {
         std::cout << "Bug in model extraction" << std::endl;
         return false;
      } else {
         if(mapIdsFacesVolum[keyU].size()!=1){
            std::cout << "Bug in model extraction : 2" << std::endl;
         } else {

            mesh3d->GetCellPoints(mapIdsFacesVolum[keyU].at(0), npts, pts);
            tetSelected.push_back(mapIdsFacesVolum[keyU].at(0)); 
            
            //Get the tetrahedron corresponding
            if(npts!=4){
               std::cout << "Bug in model extraction : 3" << std::endl;
            } else {
               size_t numFaceFound=0;
               bool foundFace = false;
               while(!foundFace){   
                  std::vector<size_t> testVector;
                  testVector.push_back(pts[idTetFaces[3*numFaceFound]]);
                  testVector.push_back(pts[idTetFaces[3*numFaceFound+1]]);
                  testVector.push_back(pts[idTetFaces[3*numFaceFound+2]]);
                  std::sort(testVector.begin(),testVector.end());
                  if( (testVector.at(0) == faceTranslated.at(0)) && (testVector.at(1) == faceTranslated.at(1))
                        && (testVector.at(2) == faceTranslated.at(2)) ){
                     foundFace=true;
                  } else {
                     numFaceFound++;
                  }
               }
               std::pair<size_t,size_t> keyPair = std::make_pair(mapIdsFacesVolum[keyU].at(0),numFaceFound);
               assoFaceSFaceV.push_back(keyPair);
            }
         }
      }
   }

   //Tetrahedron map points
   std::vector<std::vector<size_t>> vectorTetId;
   std::map<size_t,bool> pointsTet;

   for(size_t k=0; k<assoFaceSFaceV.size(); k++){
     
      size_t numTet = assoFaceSFaceV.at(k).first;
      size_t numExtFace = assoFaceSFaceV.at(k).second;

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

      std::vector<size_t> vectorInter;
      vectorInter.push_back(idF1[0]);
      pointsTet[idF1[0]]=true;
      pointsTet[idF1[1]]=true;
      pointsTet[idF1[2]]=true;
      pointsTet[pointIDout]=true;
   
      if(MeshHelpers::checkDirectVectorOrientation(v1,v2,v3)){
      
         //OK
         vectorInter.push_back(idF1[1]);
         vectorInter.push_back(idF1[2]);
         vectorInter.push_back(pointIDout);

      } else {
      
         //Reorder
      
         vectorInter.push_back(idF1[2]);
         vectorInter.push_back(idF1[1]);
         vectorInter.push_back(pointIDout);

      }

      vectorTetId.push_back(vectorInter);

   }

   size_t numPointsTet=0;
   std::vector<size_t> vectorTetIDF;
   for(std::map<size_t,bool>::iterator it = pointsTet.begin(); it != pointsTet.end(); ++it) {
       numPointsTet++;
       vectorTetIDF.push_back(it->first);
   }

   //Create the point matrix
   cv::Mat tetrahedronMat = cv::Mat::zeros(numPointsTet,3,CV_32FC1);
   for(size_t k=0;k<numPointsTet;k++){
      double pReal[3];
      mesh3d->GetPoints()->GetPoint(vectorTetIDF.at(k),pReal);
      for(size_t i=0; i<3;i++){
         tetrahedronMat.at<float>(k,i)=float(pReal[i]);
      }
   }

   //
   for(size_t k=0; k<vectorTetId.size(); k++){
      std::vector<size_t> vectorInter;
      for(size_t i=0; i<4; i++){
         auto idFoundInter = std::find(vectorTetIDF.begin(),vectorTetIDF.end(),vectorTetId.at(k).at(i));
         size_t idFound = idFoundInter - vectorTetIDF.begin();
         vectorInter.push_back(idFound);
      }
      tetIdVector.push_back(vectorInter);
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
