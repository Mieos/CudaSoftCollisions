#include "csc/MeshStructureCollider.hpp"
#include <iostream>
#include <math.h> 

#include "kernels/csc_core.h"

//Init all the buffers
MeshStructureCollider::MeshStructureCollider(const cv::Mat & dataMesh, const std::vector<std::vector<size_t> > & tetIdVector, const std::vector<size_t> & associationVectorU) : initialized(false), numPoints(0), verbose(true){

   //Cuda error
   bool errorDetected = false;

   std::vector<float> array;

   if(dataMesh.type()==CV_32FC1){

      this->numPoints=dataMesh.rows;
      this->useSubdivision=false;
      this->numberOfSpatialSubdivision=10;
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
         if(!MeshStructureCollider::cudaAllocation((void **) &(this->data_d), size_data, "Data buffer")){
            errorDetected=true;
         } else {
            if(!MeshStructureCollider::copyTodevice(this->data_d,this->dataArrayBuff,size_data,"Data buffer")){
               errorDetected=true;
            }
         }

         //Copy the index of the tetrahedrons
         this->numTets=tetIdVector.size();
         size_t* tetVectorPointer = new size_t[4*tetIdVector.size()];

         for(size_t k=0; k<tetIdVector.size(); k++){

            if(tetIdVector.at(k).size()!=4){
               std::cout << "Issue with association vector in collider initialisation, abort.."<< std::endl;
               errorDetected=true;
               break;
            }
            tetVectorPointer[4*k] = tetIdVector.at(k).at(0);
            tetVectorPointer[4*k+1] = tetIdVector.at(k).at(1);
            tetVectorPointer[4*k+2] = tetIdVector.at(k).at(2);
            tetVectorPointer[4*k+3] = tetIdVector.at(k).at(3);

         }

         //Copy index to GPUs
         size_t size_tetVector = 4*this->numTets*sizeof(size_t);
         if(!MeshStructureCollider::cudaAllocation((void **) &(this->tetId_d), size_tetVector, "Structure")){
            errorDetected=true;
         } else {
            if(!MeshStructureCollider::copyTodevice(this->tetId_d,tetVectorPointer,size_tetVector,"Structure")){
               errorDetected=true;
            }
         }

         //Copy the association vector
         this->associationVector = associationVectorU;

         //Create sphere buffer: (x y z radius)
         size_t size_sphereBuffer = 4*this->numTets*sizeof(float);
         if(!MeshStructureCollider::cudaAllocation((void **) &(this->sphereBuf_d), size_sphereBuffer, "Spheres")){
            errorDetected=true;
         }

         //Create normal buffer for tetrahedron (n1 n2 n3 n4 P0)
         size_t size_normalBuffer = 4*6*this->numTets*sizeof(float);
         if(!MeshStructureCollider::cudaAllocation((void **) &(this->normalBuf_d), size_normalBuffer, "Normals")){
            errorDetected=true;
         }

         //Create collision vector (CPU and GPU)
         this->collideVectorArray = new bool[this->numTets];
         size_t size_collideVectorArray = this->numTets*sizeof(bool); 
         if(!MeshStructureCollider::cudaAllocation((void **) &(this->collideVectorArray_d), 
                  size_collideVectorArray, "Collision vector")){
            errorDetected=true;
         }

         //Create inversion vector (CPU and GPU)
         this->inversionTetVectorArray = new bool[this->numTets];
         size_t size_inversionTetVectorArray = this->numTets*sizeof(bool); 
         if(!MeshStructureCollider::cudaAllocation((void **) &(this->inversionTetVectorArray_d), 
                  size_inversionTetVectorArray, "Inversion vector")){
            errorDetected=true;
         }

         //Subdivision buffers
         if(this->useSubdivision){
            size_t size_subdividedCollisionVector = this->numberSub*this->numTets*sizeof(bool);
            if(!MeshStructureCollider::cudaAllocation((void **) &(this->subdividedCollisionVector_d), 
                     size_subdividedCollisionVector, "Loop subdivision vector")){
               errorDetected=true;
            }
         }

         //Spatial subdivision
         if(this->useSpatialSubdivision){

            //Limits
            this->subdivisionXYZArray=new float[3*(this->numberOfSpatialSubdivision+1)];
            size_t size_subarrayBorder = 3*(this->numberOfSpatialSubdivision+1)*sizeof(float);
            if(!MeshStructureCollider::cudaAllocation((void **) &(this->subdivisionXYZArray_d), 
                     size_subarrayBorder, "Loop subdivision vector")){
               errorDetected=true;
            }

            //State of each tetrahedron
            size_t size_spatialSub = this->numberOfSpatialSubdivision*this->numberOfSpatialSubdivision*this->numberOfSpatialSubdivision*this->numTets*sizeof(bool);
            if(!MeshStructureCollider::cudaAllocation((void **) &(this->spatialSub_d), 
                     size_spatialSub, "Loop subdivision vector")){
               errorDetected=true;
            }

         }

         //TODO
         size_t size_movementsArray = 3*this->numTets*sizeof(float);
         this->movementsArray = new float[3*this->numTets];
         if(!MeshStructureCollider::cudaAllocation((void **) &(this->movementsArray_d), 
                     size_movementsArray, "Movement array")){
            errorDetected=true;
         }

         //end TODO

         //Update state
         this->initialized=!errorDetected;

         //Delete intermediate vector
         delete[] tetVectorPointer;

      } else {

         std::cout << "Wrong number of collumns" << std::endl;

      }

   } else { 
      std::cout << "Type not taken into account in MeshStructureCollider, please convert to CV_32FC1" << std::endl;
   }

}

//Cuda Allocation
bool MeshStructureCollider::cudaAllocation(void ** pointer_d, size_t sizeUsed, std::string errorName){

   //Error
   cudaError_t error;

   error=cudaMalloc(pointer_d, sizeUsed);
   if (error != cudaSuccess) {
      std::string errorString = "Issue during buffer initialisation : " + errorName + "\n";
      std::cerr << errorString;
      fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(error));
      return false;
   } else {
      return true;
   }

}

//Cuda CopyToDevice
bool MeshStructureCollider::copyTodevice(void* pointer_gpu, void* pointer_cpu, size_t sizeUsed , std::string errorName){

   //Error
   cudaError_t error;

   error=cudaMemcpy(pointer_gpu, pointer_cpu, sizeUsed, cudaMemcpyHostToDevice);
   if (error != cudaSuccess) {
      std::string errorString = "Issue loading CPU to GPU : " + errorName + "\n";
      std::cerr << errorString;
      fprintf(stderr, "cudaMemcpy failed: %s\n", cudaGetErrorString(error));
      return false;
   } else {
      return true;
   }

}

//Cuda CopyToHost
bool MeshStructureCollider::copyTohost(void* pointer_cpu, void* pointer_gpu, size_t sizeUsed , std::string errorName){

   //Error
   cudaError_t error;

   error=cudaMemcpy(pointer_cpu, pointer_gpu, sizeUsed, cudaMemcpyDeviceToHost);
   if (error != cudaSuccess) {
      std::string errorString = "Issue loading GPU to CPU : " + errorName + "\n";
      std::cerr << errorString;
      fprintf(stderr, "cudaMemcpy failed: %s\n", cudaGetErrorString(error));
      return false;
   } else {
      return true;
   }

}


bool MeshStructureCollider::checkGPUerrors(std::string errorName){
   if ( cudaSuccess != cudaGetLastError() ){ 
      std::cerr << "Error Kernel : " << errorName << std::endl;
      return false;
   }
   cudaDeviceSynchronize();
   if ( cudaSuccess != cudaGetLastError() ){ 
      std::cerr << "Error Kernel : " << errorName << " (sync)" << std::endl;
      return false;
   }
   return true;
}

//Destroy all the buffers
MeshStructureCollider::~MeshStructureCollider(){

   //CPU arrays
   delete[] this->dataArrayBuff;
   delete[] this->collideVectorArray;
   delete[] this->inversionTetVectorArray;

   //GPU arrays
   cudaFree(this->data_d);
   cudaFree(this->tetId_d);
   cudaFree(this->sphereBuf_d);
   cudaFree(this->normalBuf_d);
   cudaFree(this->inversionTetVectorArray_d);
   cudaFree(this->collideVectorArray_d);

   //Optional arrays (CPU and GPU)
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

   //Copy to GPU
   size_t size_data = 3*this->numPoints*sizeof(float); 
   if(!MeshStructureCollider::copyTodevice(this->data_d,this->dataArrayBuff,size_data,"Data Update")){
      return false;
   }

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
      if(!MeshStructureCollider::copyTodevice(this->subdivisionXYZArray_d, this->subdivisionXYZArray,
               size_subdivisionArray, "Spatial Sub Update")){
         return false;
      }

      //Update space subdivision tet
      updateSpatialSub<<<dimGrid, dimBlock>>>(this->data_d, this->tetId_d, this->numTets, this->subdivisionXYZArray_d, this->numberOfSpatialSubdivision, this->spatialSub_d);
      if(!MeshStructureCollider::checkGPUerrors("Update Spatial Sub")){
         return false;
      }

   }

   //Update sphere
   updateCenterSphere<<<dimGrid, dimBlock>>>(this->data_d, this->tetId_d, this->sphereBuf_d, this->numTets);
   if(!MeshStructureCollider::checkGPUerrors("Update Spheres")){
      return false;
   }

   //Update Normals Buff
   updatePlanesTet<<<dimGrid, dimBlock>>>(this->data_d, this->tetId_d, this->normalBuf_d, this->numTets);
   if(!MeshStructureCollider::checkGPUerrors("Update Faces buffers")){
      return false;
   }

   return true;

}

//Collision function

bool MeshStructureCollider::collide(){

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
   checkTetOrientations<<<dimGrid, dimBlock>>>(this->data_d, this->tetId_d, this->numTets, this->inversionTetVectorArray_d);
   if(!MeshStructureCollider::checkGPUerrors("Orientation Check")){
      return false;
   }

   //Check collision
   if(this->useSubdivision){
      size_t size_loop = size_t(float(this->numTets)/float(this->numberSub));
      checkForIntersectionV1<<<dimGrid, dimBlock>>>(this->data_d, this->tetId_d, this->sphereBuf_d, this->normalBuf_d, this->numTets, size_loop, this->subdividedCollisionVector_d); 
   } else {
      //TODO new
      //checkForIntersectionV0<<<dimGrid, dimBlock>>>(this->data_d, this->tetId_d, this->sphereBuf_d, this->normalBuf_d, this->numTets, this->collideVectorArray_d); 
      checkForIntersectionV0_withmovements<<<dimGrid, dimBlock>>>(this->data_d, this->tetId_d, this->sphereBuf_d, this->normalBuf_d, this->numTets, this->collideVectorArray_d, this->inversionTetVectorArray_d, this->movementsArray_d); 
   }
   if(!MeshStructureCollider::checkGPUerrors("Intersection checks")){
      return false;
   }

   //Reduce if needed (if subdivision is used)
   if(this->useSubdivision){

      dim3 dimGridReduction(sizeGridx,1);
      reduceIntersectionVector<<<dimGridReduction, dimBlock>>>(this->numTets, this->numberSub, this->subdividedCollisionVector_d, this->collideVectorArray_d);
      if(!MeshStructureCollider::checkGPUerrors("Reduction")){
         return false;
      }

   }

   return true;

}

//Colliding function
bool MeshStructureCollider::collide(std::vector<bool> & collisionList){

   //Collide
   this->collide();

   //Get results
   size_t sizeResult = this->numTets * sizeof(bool);
   if(!MeshStructureCollider::copyTohost(this->inversionTetVectorArray, this->inversionTetVectorArray_d, 
            sizeResult, "Inversion vector Copy")){
      return false;
   }
   if(!MeshStructureCollider::copyTohost(this->collideVectorArray, this->collideVectorArray_d, 
            sizeResult, "Collision vector Copy")){
      return false;
   }

   for(size_t k=0; k<this->numTets;k++){
      collisionList.at(k)= (inversionTetVectorArray[k] || collideVectorArray[k]);
   }

   return true; 

}

//Colliding function
bool MeshStructureCollider::collideAndGetMovements(std::vector<bool> & collisionList, std::vector<std::vector<float> > & movVect){

   //Collide
   this->collide();

   //Get results
   size_t sizeResult = this->numTets * sizeof(bool);
   size_t sizeMov = 3*this->numTets*sizeof(float);

   if(!MeshStructureCollider::copyTohost(this->inversionTetVectorArray, this->inversionTetVectorArray_d, 
            sizeResult, "Inversion vector Copy")){
      return false;
   }
   if(!MeshStructureCollider::copyTohost(this->collideVectorArray, this->collideVectorArray_d, 
            sizeResult, "Collision vector Copy")){
      return false;
   }
   if(!MeshStructureCollider::copyTohost(this->movementsArray, this->movementsArray_d, 
            sizeMov, "Movement vector Copy")){
      return false;
   }

   for(size_t k=0; k<this->numTets;k++){
      collisionList.at(k)= (inversionTetVectorArray[k] || collideVectorArray[k]);
      if(collisionList.at(k)){
         movVect.at(k).clear();
         movVect.at(k).push_back(movementsArray[3*k]);
         movVect.at(k).push_back(movementsArray[3*k+1]);
         movVect.at(k).push_back(movementsArray[3*k+2]);
      }
   }

   return true; 

}
