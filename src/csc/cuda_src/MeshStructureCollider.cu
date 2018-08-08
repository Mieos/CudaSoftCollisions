#include "csc/MeshStructureCollider.hpp"
#include <iostream>
#include <math.h> 

__device__ bool checkSphereIntersection(float* centerSphereID, size_t* s1, size_t* s2){


   //if((centerSphereID[4*(*s1)+3]==0)||(centerSphereID[4*(*s2)+3]==0)){
   //printf("WTF\n");
   //}

   float distCenter = 
      (centerSphereID[4*(*s1)] - centerSphereID[4*(*s2)])*
      (centerSphereID[4*(*s1)] - centerSphereID[4*(*s2)])+
      (centerSphereID[4*(*s1)+1] - centerSphereID[4*(*s2)+1])*
      (centerSphereID[4*(*s1)+1] - centerSphereID[4*(*s2)+1])+
      (centerSphereID[4*(*s1)+2] - centerSphereID[4*(*s2)+2])*
      (centerSphereID[4*(*s1)+2] - centerSphereID[4*(*s2)+2]);


   distCenter=sqrt(distCenter);
   float addRadius = centerSphereID[4*(*s1)+3] + centerSphereID[4*(*s2)+3]; 

   if(distCenter<addRadius){
      return true;
   } else {
      return false;
   }

} 

__device__ bool checkOrientation(float* Ptest, float* Nm, float* Pu){

   float vt2[3];
   for(size_t k=0; k<3; k++){
      vt2[k]=Ptest[k]-Pu[k];
   }

   float dotP = Nm[0]*vt2[0] + Nm[1]*vt2[1] + Nm[2]*vt2[2];

   //printf("%f\n",dotP);

   return (dotP>=0);

}

__device__ bool checkTetraIntersection(float*  dataPointsD, size_t* idArrayD, float* normalBuf, size_t* s1, size_t* s2){

   size_t id1T1 = idArrayD[4*(*s1)];
   size_t id2T1 = idArrayD[4*(*s1)+1];
   size_t id3T1 = idArrayD[4*(*s1)+2];
   size_t id4T1 = idArrayD[4*(*s1)+3];

   size_t numTet = (*s1);

   /*
   printf("OK2:%u\n",numTet);

   printf("px buf:%f\n",normalBuf[4*6*numTet+3]);
   printf("py buf:%f\n",normalBuf[4*6*numTet+4]);
   printf("pz buf:%f\n",normalBuf[4*6*numTet+5]);
   
   printf("nx buf:%f\n",normalBuf[4*6*numTet]);
   printf("ny buf:%f\n",normalBuf[4*6*numTet+1]);
   printf("nz buf:%f\n",normalBuf[4*6*numTet+2]);
   */

   //DEBUG test only one tet
   for(size_t k=0; k<1; k++){

      float normalU[3];
      float pU[3];
      normalU[0] = normalBuf[4*6*numTet+k*6];
      normalU[1] = normalBuf[4*6*numTet+k*6+1];
      normalU[2] = normalBuf[4*6*numTet+k*6+2];
      pU[0] = normalBuf[4*6*numTet+k*6+3];
      pU[1] = normalBuf[4*6*numTet+k*6+4];
      pU[2] = normalBuf[4*6*numTet+k*6+5];

      //DEBUG
      float pDebug[3];
      pDebug[0] = dataPointsD[3*id4T1];
      pDebug[1] = dataPointsD[3*id4T1+1];
      pDebug[2] = dataPointsD[3*id4T1+2];

      //BUG FIXME
      if(checkOrientation(pDebug,normalU,pU)){
         return true;
      }

   }

   return false;
}

__global__ void updateCenterSphere(float*  dataPointsD, size_t* idArrayD, float* sphereBuffer, size_t numberTets){

   size_t numTet = blockIdx.x*blockDim.x*blockDim.y +blockDim.x*threadIdx.y+threadIdx.x;

   size_t id1 = idArrayD[4*numTet];
   size_t id2 = idArrayD[4*numTet+1];
   size_t id3 = idArrayD[4*numTet+2];
   size_t id4 = idArrayD[4*numTet+3];

   if(numTet<numberTets){

      //Coordinate
      sphereBuffer[4*numTet] = (dataPointsD[3*id1] + dataPointsD[3*id2] + dataPointsD[3*id3] + dataPointsD[3*id4])/4.0;
      sphereBuffer[4*numTet+1] = (dataPointsD[3*id1+1] + dataPointsD[3*id2+1] + dataPointsD[3*id3+1] + dataPointsD[3*id4+1])/4.0;
      sphereBuffer[4*numTet+2] = (dataPointsD[3*id1+2] + dataPointsD[3*id2+2] + dataPointsD[3*id3+2] + dataPointsD[3*id4+2])/4.0;

      //Radius
      float radius = 
         (sphereBuffer[4*numTet] - dataPointsD[3*id1])*(sphereBuffer[4*numTet] - dataPointsD[3*id1])+
         (sphereBuffer[4*numTet+1] - dataPointsD[3*id1+1])*(sphereBuffer[4*numTet+1] - dataPointsD[3*id1+1])+
         (sphereBuffer[4*numTet+2] - dataPointsD[3*id1+2])*(sphereBuffer[4*numTet+2] - dataPointsD[3*id1+2]);
      sphereBuffer[4*numTet+3] = sqrt(radius);

      /*
      if(numTet==1){
         printf("radius = %f\n",radius);
         printf("data debug x = %f\n",dataPointsD[3*id1]);
         printf("data debug y = %f\n",dataPointsD[3*id1+1]);
         printf("data debug z = %f\n",dataPointsD[3*id1+2]);
      }
      */

   } else {
      //printf ("DEBUG\n");
   }

}

__global__ void updatePlanesTet(float*  dataPointsD, size_t* idArrayD, float* normalBuf, size_t numberTets){

   size_t numTet = blockIdx.x*blockDim.x*blockDim.y +blockDim.x*threadIdx.y+threadIdx.x;

   size_t idTetConsidered[4];
   idTetConsidered[0] = idArrayD[4*numTet];
   idTetConsidered[1] = idArrayD[4*numTet+1];
   idTetConsidered[2] = idArrayD[4*numTet+2];
   idTetConsidered[3] = idArrayD[4*numTet+3];

   size_t idOrientedTetFaces[] = {
      0,1,2,
      0,2,3,
      0,3,1,
      2,1,3
   };

   if(numTet<numberTets){

      for(size_t k=0; k<4; k++){

      
         size_t id1,id2,id3;
         id1=idTetConsidered[idOrientedTetFaces[3*k]];
         id2=idTetConsidered[idOrientedTetFaces[3*k+1]];
         id3=idTetConsidered[idOrientedTetFaces[3*k+2]];

         //Point of the face (oriented)
         float P1[3];
         float P2[3];
         float P3[3];
         for(size_t i=0; i<3; i++){
            P1[i] = dataPointsD[3*id1+i];
            P2[i] = dataPointsD[3*id2+i];
            P3[i] = dataPointsD[3*id3+i];
         }

         //Vectors
         float v1[3], v2[3];
         for(size_t i=0; i<3; i++){
            v1[i]=P2[i]-P1[i];
            v2[i]=P3[i]-P1[i];
         }

         //Vect product
         float v1Vv2[3];
         v1Vv2[0] = v1[1]*v2[2] - v1[2]*v2[1];
         v1Vv2[1] = v1[2]*v2[0] - v1[0]*v2[2];
         v1Vv2[2] = v1[0]*v2[1] - v1[1]*v2[0];

         //Update buffer
         float normVP = v1Vv2[0]*v1Vv2[0] + v1Vv2[1]*v1Vv2[1] + v1Vv2[2]*v1Vv2[2];
         normVP = sqrt(normVP);

         //Update buffer
         normalBuf[4*6*numTet+k*6] = v1Vv2[0]/normVP; //Normal
         normalBuf[4*6*numTet+k*6+1] = v1Vv2[1]/normVP; 
         normalBuf[4*6*numTet+k*6+2] = v1Vv2[2]/normVP;
         normalBuf[4*6*numTet+k*6+3] = dataPointsD[3*id1]; //Point
         normalBuf[4*6*numTet+k*6+4] = dataPointsD[3*id1+1];
         normalBuf[4*6*numTet+k*6+5] = dataPointsD[3*id1+2];

         /*
         if((numTet==1)&&(k==0)){

            printf("DEBUG update px = %f\n",normalBuf[4*6*numTet+k*6+3]);
            printf("DEBUG update py = %f\n",normalBuf[4*6*numTet+k*6+4]);
            printf("DEBUG update pz = %f\n",normalBuf[4*6*numTet+k*6+5]);
 
         }
         */

      }

   } else {
      //printf ("DEBUG\n");
   }

}

__global__ void checkForIntersection(float*  dataPointsD, size_t* idArrayD, float* centerSphereB, float* normalsB, size_t numberTets, bool* intersectionVector){

   size_t numTet = blockIdx.x*blockDim.x*blockDim.y +blockDim.x*threadIdx.y+threadIdx.x;
   intersectionVector[numTet]=false;

   //DEBUG
   //if(numTet==2){

      //printf("OK1\n");

      if(numTet<numberTets){
         //First part
         for(size_t k=0; k<numTet-1; k++){
            if(checkSphereIntersection(centerSphereB,&numTet,&k)){
               //printf("DEB : %u\n", numTet);
               if(checkTetraIntersection(dataPointsD,idArrayD, normalsB, &numTet,&k)){
                  intersectionVector[numTet]=true;
                  break;
               }
            }
         }
         //Second part
         if(!intersectionVector[numTet]){
            for(size_t k=numTet+1; k<numberTets; k++){ 
               if(checkSphereIntersection(centerSphereB,&numTet,&k)){ 
                  if(checkTetraIntersection(dataPointsD,idArrayD,normalsB, &numTet,&k)){
                     intersectionVector[numTet]=true;
                  }
               }
            }
         }
      } else { 
         //printf ("DEBUG\n");
      }

   //}

}

MeshStructureCollider::MeshStructureCollider(const cv::Mat & dataMesh, const std::vector<std::vector<size_t> > & tetIdVector, const std::vector<size_t> & associationVectorU) : initialized(false), numPoints(0), verbose(true){

   std::vector<float> array;

   if(dataMesh.type()==CV_32FC1){

      //std::cout << "Type detected : CV_32FC1 : OK" << std::endl;

      this->numPoints=dataMesh.rows;

      if(dataMesh.cols==3){

         this->dataArrayBuff = new float[3*numPoints];

         //Copy the points
         for(size_t i=0; i<dataMesh.rows; i++){
            for(size_t j=0; j<3; j++){
               this->dataArrayBuff[3*i+j]=dataMesh.at<float>(i,j);   
            }
         }

         //Copy points to gpu
         size_t size_data = 3*this->numPoints*sizeof(float);
         cudaMalloc((void **) &(this->data_d), size_data);
         cudaMemcpy(this->data_d, this->dataArrayBuff, size_data, cudaMemcpyHostToDevice);

         //Copy the index of the tetrahedrons
         this->numTets=tetIdVector.size();
         size_t* tetVectorPointer = new size_t[4*tetIdVector.size()];

         for(size_t k=0; k<tetIdVector.size(); k++){

            if(tetIdVector.at(k).size()!=4){
               std::cout << "Issue with association vector in collider initialisation, abort.."<< std::endl;
               break;
            }

            tetVectorPointer[4*k] = tetIdVector.at(k).at(0);
            tetVectorPointer[4*k+1] = tetIdVector.at(k).at(1);
            tetVectorPointer[4*k+2] = tetIdVector.at(k).at(2);
            tetVectorPointer[4*k+3] = tetIdVector.at(k).at(3);

         }

         size_t size_tetVector = 4*this->numTets*sizeof(size_t);
         cudaMalloc((void **) &(this->tetId_d), size_tetVector);
         cudaMemcpy(this->tetId_d, tetVectorPointer, size_tetVector, cudaMemcpyHostToDevice);

         //Copy the association vector
         this->associationVector = associationVectorU;

         //Create sphere buffer: (x y z radius)
         size_t size_sphereBuffer = 4*this->numTets*sizeof(float);
         cudaMalloc((void **) &(this->sphereBuf_d), size_sphereBuffer);

         //Create normal buffer for tetrahedron (n1 n2 n3 n4 P0)
         size_t size_normalBuffer = 4*6*this->numTets*sizeof(float);
         cudaMalloc((void **) &(this->normalBuf_d), size_normalBuffer);

         //Create collision vector (CPU and GPU)
         this->collideVectorArray = new bool[this->numTets];
         size_t size_collideVectorArray = this->numTets*sizeof(bool); 
         cudaMalloc((void **) &(this->collideVectorArray_d), size_collideVectorArray);

         //Update state
         this->initialized=true;

         delete[] tetVectorPointer;

      } else {

         std::cout << "Wrong number of collumns" << std::endl;

      }

   } else { 
      std::cout << "Type not taken into account in MeshStructureCollider, please convert to CV_32FC1" << std::endl;
   }

}

MeshStructureCollider::~MeshStructureCollider(){

   delete[] this->dataArrayBuff;
   delete[] this->collideVectorArray;
   cudaFree(this->data_d);
   cudaFree(this->tetId_d);
   cudaFree(this->sphereBuf_d);
   cudaFree(this->normalBuf_d);

}

bool MeshStructureCollider::isProperlyInitialized(){

   return this->initialized;

}

bool MeshStructureCollider::updatePointsPositions(const cv::Mat & newPositions){

   //Not init
   if(!this->isProperlyInitialized()){
      return false;
   }

   //Bad number of cols
   if(newPositions.cols!=3){
      return false;
   }

   //Test type
   if(newPositions.type()==CV_32FC1){

      for(size_t k=0; k<this->associationVector.size(); k++){
         this->dataArrayBuff[3*k]=newPositions.at<float>(associationVector.at(k),0); 
         this->dataArrayBuff[3*k+1]=newPositions.at<float>(associationVector.at(k),1); 
         this->dataArrayBuff[3*k+2]=newPositions.at<float>(associationVector.at(k),2); 
      }

   } else if(newPositions.type()==CV_64FC1){

      for(size_t k=0; k<this->associationVector.size(); k++){
         this->dataArrayBuff[3*k]=float(newPositions.at<double>(associationVector.at(k),0)); 
         this->dataArrayBuff[3*k+1]=float(newPositions.at<double>(associationVector.at(k),1)); 
         this->dataArrayBuff[3*k+2]=float(newPositions.at<double>(associationVector.at(k),2)); 
      }

   } 

   //Copy to buffer
   size_t size_data = 3*this->numPoints*sizeof(float);
   cudaMemcpy(this->data_d, this->dataArrayBuff, size_data, cudaMemcpyHostToDevice);

   return true;

}

bool MeshStructureCollider::collide(std::vector<bool> & collisionList){

   cudaDeviceProp prop;
   cudaGetDeviceProperties(&prop,0);
   size_t sizeBlockToUse = sqrt (prop.maxThreadsDim[0]);
   size_t sizeGridx = ceil(float(this->numTets)/float(sizeBlockToUse*sizeBlockToUse));
   size_t sizeGridy = 1;

   //std::cout << "Tets = " << this->numTets << std::endl;
   //std::cout << "Block size = " << sizeBlockToUse << std::endl;
   //std::cout << "Grid size = " << sizeGridx << std::endl;

   dim3 dimGrid(sizeGridx,sizeGridy);
   dim3 dimBlock(sizeBlockToUse, sizeBlockToUse);

   //Update sphere
   updateCenterSphere<<<dimGrid, dimBlock>>>(this->data_d, this->tetId_d, this->sphereBuf_d, this->numTets);
   cudaDeviceSynchronize();

   //Update Normals Buff
   updatePlanesTet<<<dimGrid, dimBlock>>>(this->data_d, this->tetId_d, this->normalBuf_d, this->numTets);
   cudaDeviceSynchronize();

   //Check collision
   checkForIntersection<<<dimGrid, dimBlock>>>(this->data_d, this->tetId_d, this->sphereBuf_d, this->normalBuf_d, this->numTets, this->collideVectorArray_d); 
   cudaDeviceSynchronize();

   //Get results
   size_t sizeResult = sizeof(bool)*this->numTets;
   cudaMemcpy(collideVectorArray, collideVectorArray_d, sizeResult, cudaMemcpyDeviceToHost);
   for(size_t k=0; k<this->numTets;k++){
      collisionList.at(k)=collideVectorArray[k];
   }

   return false; 

}
