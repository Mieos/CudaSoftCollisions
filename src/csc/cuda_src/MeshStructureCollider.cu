#include "csc/MeshStructureCollider.hpp"
#include <iostream>
#include <math.h> 

#include "kernels/csc_core.h"

//Init all the buffers
MeshStructureCollider::MeshStructureCollider(const cv::Mat & dataMesh, const std::vector<std::vector<size_t> > & tetIdVector, const std::vector<size_t> & associationVectorU) : initialized(false), numPoints(0), verbose(true){

   //Cuda error
   cudaError_t error;

   std::vector<float> array;

   if(dataMesh.type()==CV_32FC1){

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
         error=cudaMalloc((void **) &(this->data_d), size_data);
         if (error != cudaSuccess) {
            fprintf(stderr, "Issue while initializing data buffer\n");
            fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(error));
         }
         error=cudaMemcpy(this->data_d, this->dataArrayBuff, size_data, cudaMemcpyHostToDevice);
         if (error != cudaSuccess) {
            fprintf(stderr, "Issue while loading data on GPU\n");
            fprintf(stderr, "cudaMemcpy failed: %s\n", cudaGetErrorString(error));
         }


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
         error=cudaMalloc((void **) &(this->tetId_d), size_tetVector);
         if (error != cudaSuccess) {
            fprintf(stderr, "Issue while initializing data buffer (structure)\n");
            fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(error));
         }

         error=cudaMemcpy(this->tetId_d, tetVectorPointer, size_tetVector, cudaMemcpyHostToDevice);
         if (error != cudaSuccess) {
            fprintf(stderr, "Issue while loading data (structure) on GPU\n");
            fprintf(stderr, "cudaMemcpy failed: %s\n", cudaGetErrorString(error));
         }

         //Copy the association vector
         this->associationVector = associationVectorU;

         //Create sphere buffer: (x y z radius)
         size_t size_sphereBuffer = 4*this->numTets*sizeof(float);
         error=cudaMalloc((void **) &(this->sphereBuf_d), size_sphereBuffer);
         if (error != cudaSuccess) {
            fprintf(stderr, "Issue while initializing buffer of spheres\n");
            fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(error));
         }


         //Create normal buffer for tetrahedron (n1 n2 n3 n4 P0)
         size_t size_normalBuffer = 4*6*this->numTets*sizeof(float);
         error=cudaMalloc((void **) &(this->normalBuf_d), size_normalBuffer);
         if (error != cudaSuccess) {
            fprintf(stderr, "Issue while initializing normal buffer\n");
            fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(error));
         }

         //Create collision vector (CPU and GPU)
         this->collideVectorArray = new bool[this->numTets];
         size_t size_collideVectorArray = this->numTets*sizeof(bool); 
         error=cudaMalloc((void **) &(this->collideVectorArray_d), size_collideVectorArray);
         if (error != cudaSuccess) {
            fprintf(stderr, "Issue while creating output vector\n");
            fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(error));
         }

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

//Destroy all the buffers
MeshStructureCollider::~MeshStructureCollider(){

   delete[] this->dataArrayBuff;
   delete[] this->collideVectorArray;
   cudaFree(this->data_d);
   cudaFree(this->tetId_d);
   cudaFree(this->sphereBuf_d);
   cudaFree(this->normalBuf_d);

}

//Check if all the buffers has been properly initialized (Not used: TODO)
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

//Colliding function
bool MeshStructureCollider::collide(std::vector<bool> & collisionList){

   cudaError_t error;

   //Get the max of threads usable
   cudaDeviceProp prop;
   cudaGetDeviceProperties(&prop,0);
   size_t sizeBlockToUse = sqrt (prop.maxThreadsDim[0]);
   
   //Reduce the number of thread in one direction in order to avoid ressources exhaustment
   size_t dimXReducedBlock = size_t(0.75 * float(sizeBlockToUse)); 
   
   //Compute size of the grid needed
   size_t sizeGridx = ceil(float(this->numTets)/float(dimXReducedBlock*sizeBlockToUse));
   size_t sizeGridy = 1;

   dim3 dimGrid(sizeGridx,sizeGridy);
   dim3 dimBlock(dimXReducedBlock, sizeBlockToUse);

   //Update sphere
   updateCenterSphere<<<dimGrid, dimBlock>>>(this->data_d, this->tetId_d, this->sphereBuf_d, this->numTets);
   if ( cudaSuccess != cudaGetLastError() ){ 
      fprintf(stderr, "Error updating spheres\n");
   }
   cudaDeviceSynchronize();
   if ( cudaSuccess != cudaGetLastError() ){ 
      fprintf(stderr, "Error updating spheres (synchro)\n");
   }

   //Update Normals Buff
   updatePlanesTet<<<dimGrid, dimBlock>>>(this->data_d, this->tetId_d, this->normalBuf_d, this->numTets);
   if ( cudaSuccess != cudaGetLastError() ){ 
      fprintf(stderr, "Error updating face buffers\n");
   }
   cudaDeviceSynchronize();
   if ( cudaSuccess != cudaGetLastError() ){ 
      fprintf(stderr, "Error updating face buffers(synchro)\n");
   }

   //Check collision
   checkForIntersection<<<dimGrid, dimBlock>>>(this->data_d, this->tetId_d, this->sphereBuf_d, this->normalBuf_d, this->numTets, this->collideVectorArray_d); 
   error=cudaGetLastError();
   if ( cudaSuccess != error ){ 
      fprintf(stderr, "Error collision\n");
      fprintf(stderr, "Kernel: %s\n", cudaGetErrorString(error));
   }
   error = cudaDeviceSynchronize();
   if ( cudaSuccess != error ){ 
      fprintf(stderr, "Error collision(synchro)\n");
      fprintf(stderr, "Synchro: %s\n", cudaGetErrorString(error));
   }

   //Get results
   size_t sizeResult = this->numTets * sizeof(bool);

   error = cudaMemcpy(collideVectorArray, this->collideVectorArray_d, sizeResult, cudaMemcpyDeviceToHost);
   if (error != cudaSuccess) {
      fprintf(stderr, "Issue while saving the collision vector\n");
      fprintf(stderr, "cudaMemcpy failed: %s\n", cudaGetErrorString(error));
   }

   for(size_t k=0; k<this->numTets;k++){
      collisionList.at(k)=collideVectorArray[k];
   }

   return false; 

}
