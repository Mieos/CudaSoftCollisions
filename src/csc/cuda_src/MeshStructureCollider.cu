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
      this->useSubdivision=false;
      this->useSpatialSubdivision=false; //NOT WORKING FIXME

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

         //Subdivision buffers
         if(this->useSubdivision){
            this->numberSub=10;
            size_t size_subdividedCollisionVector = this->numberSub*this->numTets*sizeof(bool);
            error=cudaMalloc((void **) &(this->subdividedCollisionVector_d), size_subdividedCollisionVector);
            if (error != cudaSuccess) {
               fprintf(stderr, "Issue while creating subdivision vector\n");
               fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(error));
            }
         }

         //Create collision vector (CPU and GPU)
         this->collideVectorArray = new bool[this->numTets];
         size_t size_collideVectorArray = this->numTets*sizeof(bool); 
         error=cudaMalloc((void **) &(this->collideVectorArray_d), size_collideVectorArray);
         if (error != cudaSuccess) {
            fprintf(stderr, "Issue while creating output vector\n");
            fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(error));
         }

         //Spatial subdivision
         if(this->useSpatialSubdivision){

            //Limits
            this->numberOfSpatialSubdivision=10;
            this->subdivisionXYZArray=new float[3*(this->numberOfSpatialSubdivision+1)];
            size_t size_subarrayBorder = 3*(this->numberOfSpatialSubdivision+1)*sizeof(float);
            error=cudaMalloc((void **) &(this->subdivisionXYZArray_d), size_subarrayBorder);
            if (error != cudaSuccess) {
               fprintf(stderr, "Issue while creating spatial subdivision vector\n");
               fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(error));
            }

            //State of each tetrahedron
            size_t size_spatialSub = this->numberOfSpatialSubdivision*this->numberOfSpatialSubdivision*this->numberOfSpatialSubdivision*this->numTets*sizeof(bool);
            error=cudaMalloc((void **) &(this->spatialSub_d), size_spatialSub);
            if (error != cudaSuccess) {
               fprintf(stderr, "Issue while creating spatial subdivision vector (2)\n");
               fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(error));
            }

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

   if(this->useSubdivision){
      cudaFree(this->subdividedCollisionVector_d);
   }

   if(this->useSpatialSubdivision){
      delete[] this->subdivisionXYZArray;
      cudaFree(this->spatialSub_d);
      cudaFree(this->subdivisionXYZArray_d);
   }

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

   //Compute Bounding boxes of the global space
   float maxX, maxY, maxZ, minX, minY, minZ;
   maxX = FLT_MIN;
   maxY = FLT_MIN;
   maxZ = FLT_MIN;
   minX = FLT_MAX;
   minY = FLT_MAX;
   minZ = FLT_MAX;

   //Test type
   if(newPositions.type()==CV_32FC1){

      for(size_t k=0; k<this->associationVector.size(); k++){

         this->dataArrayBuff[3*k]=newPositions.at<float>(associationVector.at(k),0); 
         this->dataArrayBuff[3*k+1]=newPositions.at<float>(associationVector.at(k),1); 
         this->dataArrayBuff[3*k+2]=newPositions.at<float>(associationVector.at(k),2);

         //Updata BB values
         if(maxX<this->dataArrayBuff[3*k]){
            maxX=this->dataArrayBuff[3*k];
         }
         if(maxY<this->dataArrayBuff[3*k+1]){
            maxY=this->dataArrayBuff[3*k+1];
         }
         if(maxZ<this->dataArrayBuff[3*k+2]){
            maxZ=this->dataArrayBuff[3*k+2];
         }
         if(minX>this->dataArrayBuff[3*k]){
            minX=this->dataArrayBuff[3*k];
         }
         if(minY>this->dataArrayBuff[3*k+1]){
            minY=this->dataArrayBuff[3*k+1];
         }
         if(minZ>this->dataArrayBuff[3*k+2]){
            minZ=this->dataArrayBuff[3*k+2];
         }

      }

   } else if(newPositions.type()==CV_64FC1){

      for(size_t k=0; k<this->associationVector.size(); k++){

         this->dataArrayBuff[3*k]=float(newPositions.at<double>(associationVector.at(k),0)); 
         this->dataArrayBuff[3*k+1]=float(newPositions.at<double>(associationVector.at(k),1)); 
         this->dataArrayBuff[3*k+2]=float(newPositions.at<double>(associationVector.at(k),2)); 

         //Updata BB values
         if(maxX<this->dataArrayBuff[3*k]){
            maxX=this->dataArrayBuff[3*k];
         }
         if(maxY<this->dataArrayBuff[3*k+1]){
            maxY=this->dataArrayBuff[3*k+1];
         }
         if(maxZ<this->dataArrayBuff[3*k+2]){
            maxZ=this->dataArrayBuff[3*k+2];
         }
         if(minX>this->dataArrayBuff[3*k]){
            minX=this->dataArrayBuff[3*k];
         }
         if(minY>this->dataArrayBuff[3*k+1]){
            minY=this->dataArrayBuff[3*k+1];
         }
         if(minZ>this->dataArrayBuff[3*k+2]){
            minZ=this->dataArrayBuff[3*k+2];
         }

      }

   }

   size_t size_data = 3*this->numPoints*sizeof(float);
   cudaMemcpy(this->data_d, this->dataArrayBuff, size_data, cudaMemcpyHostToDevice);

   if(this->useSpatialSubdivision){

      float sizeXBB = (maxX-minX)/float(this->numberOfSpatialSubdivision);
      float sizeYBB = (maxY-minY)/float(this->numberOfSpatialSubdivision);
      float sizeZBB = (maxZ-minZ)/float(this->numberOfSpatialSubdivision);

      for(size_t k=0; k<this->numberOfSpatialSubdivision+1; k++){
         this->subdivisionXYZArray[k]=minX+k*sizeXBB;
         this->subdivisionXYZArray[this->numberOfSpatialSubdivision+k]=minY+k*sizeYBB;
         this->subdivisionXYZArray[2*this->numberOfSpatialSubdivision+k]=minZ+k*sizeZBB;
      }

      //Copy to buffer
      size_t size_subdivisionArray=3*this->numberOfSpatialSubdivision*sizeof(float);
      cudaMemcpy(this->subdivisionXYZArray_d, this->subdivisionXYZArray, size_subdivisionArray, cudaMemcpyHostToDevice);

      //Update space subdivision tet
      updateSpatialSub<<<dimGrid, dimBlock>>>(this->data_d, this->tetId_d, this->numTets, this->subdivisionXYZArray_d, this->numberOfSpatialSubdivision, this->spatialSub_d);

      if ( cudaSuccess != cudaGetLastError() ){ 
         fprintf(stderr, "Error updating subdivision array\n");
      }
      cudaDeviceSynchronize();
      if ( cudaSuccess != cudaGetLastError() ){ 
         fprintf(stderr, "Error updating subdivision array (synchro)\n");
      }

   }

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
   size_t sizeGridy;
   if(this->useSubdivision){
      sizeGridy = this->numberSub;
   } else {
      sizeGridy = 1;
   }
   dim3 dimGrid(sizeGridx,sizeGridy);
   dim3 dimBlock(dimXReducedBlock, sizeBlockToUse);

   //Check orientation (over time, tetrahedron can invert, this should not happen)
   dim3 dimGridOrientation(sizeGridx,1);
   checkTetOrientations<<<dimGrid, dimBlock>>>(this->data_d, this->tetId_d, this->numTets, this->collideVectorArray_d);
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


   //Check collision
   if(this->useSubdivision){
      size_t size_loop = size_t(float(this->numTets)/float(this->numberSub));
      checkForIntersectionV1<<<dimGrid, dimBlock>>>(this->data_d, this->tetId_d, this->sphereBuf_d, this->normalBuf_d, this->numTets, size_loop, this->subdividedCollisionVector_d); 
   } else {
      checkForIntersectionV0<<<dimGrid, dimBlock>>>(this->data_d, this->tetId_d, this->sphereBuf_d, this->normalBuf_d, this->numTets, this->collideVectorArray_d); 
   }

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

   //Reduce if needed (if subdivision is used)
   if(this->useSubdivision){

      dim3 dimGridReduction(sizeGridx,1);
      reduceIntersectionVector<<<dimGridReduction, dimBlock>>>(this->numTets, this->numberSub, this->subdividedCollisionVector_d, this->collideVectorArray_d);
      error=cudaGetLastError();
      if ( cudaSuccess != error ){ 
         fprintf(stderr, "Error reduction\n");
         fprintf(stderr, "Kernel: %s\n", cudaGetErrorString(error));
      } 
      error = cudaDeviceSynchronize();
      if ( cudaSuccess != error ){ 
         fprintf(stderr, "Error reduction(synchro)\n");
         fprintf(stderr, "Synchro: %s\n", cudaGetErrorString(error));
      }

   }

   //Get results
   size_t sizeResult = this->numTets * sizeof(bool);
   error = cudaMemcpy(this->collideVectorArray, this->collideVectorArray_d, sizeResult, cudaMemcpyDeviceToHost);
   if (error != cudaSuccess) {
      fprintf(stderr, "Issue while saving the collision vector\n");
      fprintf(stderr, "cudaMemcpy failed: %s\n", cudaGetErrorString(error));
   }

   for(size_t k=0; k<this->numTets;k++){
      collisionList.at(k)=collideVectorArray[k];
   }

   return true; 

}
